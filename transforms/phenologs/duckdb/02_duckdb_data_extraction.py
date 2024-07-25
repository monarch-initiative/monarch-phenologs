import duckdb
import os

duckdb.sql("CREATE TABLE nodes AS SELECT * FROM '../../../datasets/sources/monarch_kg/monarch-kg_nodes.tsv'") # Create the nodes table
duckdb.sql("CREATE TABLE edges AS SELECT * FROM '../../../datasets/sources/monarch_kg/monarch-kg_edges.tsv'") # Create the edges table

# nodes columns: id	category	name	description	xref	provided_by	synonym	full_name	in_taxon	in_taxon_label	symbol	deprecated	iri	same_as
# edges columns: id	original_subject	predicate	original_object	category	agent_type	aggregator_knowledge_source	knowledge_level	primary_knowledge_source	qualifiers	provided_by	has_evidence	publications	stage_qualifier	frequency_qualifier	has_count
# want everything that has a 'has_phenotype' edge, preferably from just the species of interest

if not os.path.exists('../../../datasets/intermediate/duckdb_tables/'):
    os.makedirs('../../../datasets/intermediate/duckdb_tables/')


# create genes table
duckdb.sql("CREATE TABLE genes AS "
           "select distinct nodes.id as gene_id, "
           "nodes.name as gene_name, "
           "nodes.in_taxon as gene_taxon_id, "
           "nodes.in_taxon_label as gene_taxon_label "
           "from nodes "
           "where nodes.in_taxon in ('NCBITaxon:9606','NCBITaxon:10090','NCBITaxon:8355','NCBITaxon:6239','NCBITaxon:10116','NCBITaxon:8364') "
           # "where nodes.in_taxon in ('NCBITaxon:9606','NCBITaxon:8355','NCBITaxon:8364') "
           "and nodes.category = 'biolink:Gene' ")

print("Genes table sample: ")
duckdb.sql("SELECT count(*) as gene_row_count FROM genes").show(max_width=10000, max_rows=100)

duckdb.sql("COPY genes TO '../../../datasets/intermediate/duckdb_tables/genes.csv' (HEADER true, DELIMITER '\t');")



# create phenotypes table
# NOTE: Should this be amended to include the taxon, so that we definitely say phenotype X is in species Y for the rat/mouse phenotypes?
# NOTE: This approach will only include phenotypes that are a phenotype of a gene, so it will drop out phenotypes with no known causal genes.
duckdb.sql("CREATE TABLE phenotypes AS "
           "select distinct nodes_B.id as phenotype_id, "
           "nodes_B.name as phenotype_name, "
           "nodes_A.in_taxon, " # This will allow us to assign a taxon for the ambiguous mouse/rat phenotypes   
           # "row_number() over (partition by nodes_A.in_taxon) as phenotype_row_number "
           "from nodes as nodes_A left join edges on nodes_A.id = edges.subject "
           "left join nodes as nodes_B "
           "on edges.object = nodes_B.id "
           "where nodes_A.category = 'biolink:Gene' "
           "and edges.predicate = 'biolink:has_phenotype' "
           "and nodes_B.category = 'biolink:PhenotypicFeature' and ( "     
           "regexp_replace(nodes_B.id, ':.+', ':') like 'HP:%' "
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'MP:%' " # Both rat and mouse use MP
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'ZP:%' "
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'XPO:%' "
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'WBPhenotype:%' "
           ") order by nodes_B.id")

#TODO: Note that XPO phenotypes have edges for two different taxons: NCBITaxon:8364	(Xenopus tropicalis) and NCBITaxon:8355 (Xenopus laevis).
# Keep both? Are the row counts identical for each taxon? Gene associations?

duckdb.sql("CREATE TABLE phenotypes_test AS SELECT p.*, row_number() over (partition by p.in_taxon) as phenotype_row_number "
           "from phenotypes p ")

print("Phenotypes table sample:")
duckdb.sql("SELECT count(*) as phenotype_row_count FROM phenotypes").show(max_width=10000, max_rows=100)

# Creating a test/development option where we only save 10 or 100 phenotypes per taxon,
# reducing the amount of data processing necessary in downstream steps for testing purposes.
# duckdb.sql("DELETE FROM phenotypes_test WHERE phenotype_row_number > 100")

print("Total number of phenotypes included:")
duckdb.sql("SELECT count(*) as phenotype_row_count FROM phenotypes_test").show(max_width=10000, max_rows=100)

duckdb.sql("COPY phenotypes_test TO '../../../datasets/intermediate/duckdb_tables/phenotypes.csv' (HEADER true, DELIMITER '\t');")




# create orthologs table - NOTE: this table is currently bidirectional, meaning when searching for the ortholog between gene A and gene B
# you will get one row returned when searching for gene A and gene B, but you will also get one row returned if you reverse the order (gene B and gene A).
# This allows you to get the ortholog between any gene pair without needing to duplicate the query, but you also must take care not to duplicate the query.
duckdb.sql("CREATE TABLE orthologs AS "
           "select distinct gene_a_id, gene_a_taxon_id, predicate, gene_b_id, gene_b_taxon_id, ortholog_id from ( "
           "select distinct gene_a.gene_id as gene_a_id, gene_A.gene_taxon_id as gene_a_taxon_id, predicate, gene_b.gene_id as gene_b_id, gene_b.gene_taxon_id as gene_b_taxon_id, regexp_replace(has_evidence, '.+:', '')  as ortholog_id "
           "from edges inner join genes as gene_A on edges.subject = gene_a.gene_id "
           "inner join genes as gene_B on edges.object = gene_B.gene_id "
           "where edges.predicate = 'biolink:orthologous_to' "
           "union all "
           "select distinct gene_b.gene_id as gene_a_id, gene_b.gene_taxon_id as gene_a_taxon_id, predicate, gene_a.gene_id as gene_b_id, gene_a.gene_taxon_id as gene_b_taxon_id, regexp_replace(has_evidence, '.+:', '')  as ortholog_id "
           "from edges inner join genes as gene_A on edges.subject = gene_a.gene_id "
           "inner join genes as gene_B on edges.object = gene_B.gene_id "
           "where edges.predicate = 'biolink:orthologous_to' )")


print("Orthologs human-mouse sample:")
duckdb.sql("SELECT * FROM orthologs "
           "where (gene_a_id = 'HGNC:16435' and gene_b_id = 'MGI:95886') or "
           "(gene_b_id = 'HGNC:16435' and gene_a_id = 'MGI:95886') ").show(max_width=10000, max_rows=10)

print("Orthologs table  sample:")
duckdb.sql("SELECT count(*) as ortholog_row_count FROM orthologs").show(max_width=10000, max_rows=10)
# total count: 178051
# distinct count: 178030
# bidirectional via union all: 356060

# write to file
duckdb.sql("COPY orthologs TO '../../../datasets/intermediate/duckdb_tables/orthologs.csv' (HEADER true, DELIMITER '\t');")


# Create a table containing the total number of common orthologs between species pairs:
# Note: This table should be bidirectional.

duckdb.sql("CREATE TABLE total_common_orthologs AS "
           "select gene_a_taxon_id as species_a_taxon_id, gene_b_taxon_id as species_b_taxon_id, count(distinct ortholog_id) as ortholog_count "
           "from orthologs "
           "group by gene_a_taxon_id, gene_b_taxon_id "
           "order by gene_a_taxon_id, gene_b_taxon_id ")

print("Common ortholog total counts:" )
duckdb.sql("SELECT * FROM total_common_orthologs order by species_a_taxon_id, species_b_taxon_id ").show(max_width=10000, max_rows=100)

# write to file
duckdb.sql("COPY total_common_orthologs TO '../../../datasets/intermediate/duckdb_tables/total_common_orthologs.csv' (HEADER true, DELIMITER '\t');")


# create gene to phenotype table
duckdb.sql("CREATE TABLE gene_to_phenotype AS "
            "select distinct nodes_A.id as gene_id, nodes_A.name as gene_name, "
           "nodes_A.in_taxon as gene_taxon_id, "
           "nodes_A.in_taxon_label as gene_taxon_label, "
           "edges.predicate, nodes_B.id as phenotype_id, "
           "nodes_B.name as phenotype_name, "
           "nodes_A.in_taxon as phenotype_taxon_id " # Assumption: If gene in taxon has phenotype, that phenotype occurs in the same taxon.
           "from nodes as nodes_A left join edges on nodes_A.id = edges.subject "
           "left join nodes as nodes_B "
            "on edges.object = nodes_B.id "
           "where nodes_A.in_taxon in ('NCBITaxon:9606','NCBITaxon:10090','NCBITaxon:8355','NCBITaxon:6239','NCBITaxon:10116','NCBITaxon:8364') "
            "and nodes_A.category = 'biolink:Gene' "
            "and edges.predicate = 'biolink:has_phenotype' "
            "and nodes_B.category = 'biolink:PhenotypicFeature'")


print("Gene to phenotype sample:")
duckdb.sql("SELECT * as ortholog_row_count FROM gene_to_phenotype").show(max_width=10000, max_rows=10)

# write to file
duckdb.sql("COPY gene_to_phenotype TO '../../../datasets/intermediate/duckdb_tables/gene_to_phenotype.csv' (HEADER true, DELIMITER '\t');")




# create gene to ortholog table, subset to common orthologs table for taxons of interest?


# create phenotype-ortholog file
# gene to phenotype schema: gene_id	gene_name	gene_taxon_id	gene_taxon_label	predicate	phenotype_id	phenotype_name
# ortholog schema: gene_a_id	gene_a_taxon_id	predicate	gene_b_id	gene_b_taxon_id	ortholog_id
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS "
           "SELECT distinct g2p.phenotype_id, g2p.phenotype_name, g2p.phenotype_taxon_id, orth.ortholog_id FROM gene_to_phenotype as g2p left join orthologs as orth on g2p.gene_id = orth.gene_a_id "
           "where orth.ortholog_id is not null")

# write to file
duckdb.sql("COPY phenotype_to_ortholog TO '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv' (HEADER true, DELIMITER '\t');")

print("Phenotype to ortholog row count:")
duckdb.sql("SELECT count(*) as phenotype_ortholog_row_count FROM phenotype_to_ortholog").show(max_width=10000, max_rows=10)


# total count: 1999745
# distinct count: 447587 -> after adding phenotype taxon ID: 448251 -> Without the phenotype taxon, the overlapping mouse/rat phenotypes collapse.


# generate phenotype-ortholog counts
'''
Here we need to put together a simple count of distinct orthologs per phenotype
Can utilize the phenotype_to_ortholog table.
Columns: phenotype_id	phenotype_name	phenotype_taxon_id	ortholog_id
'''

duckdb.sql("CREATE TABLE phenotype_ortholog_counts AS "
           "SELECT DISTINCT phenotype_id, phenotype_taxon_id, count(distinct ortholog_id) as ortholog_count "
           "FROM phenotype_to_ortholog "
           "GROUP BY phenotype_id, phenotype_taxon_id ")
           # "where phenotype_a.phenotype_id in ('MP:0001399','MP:0001025','MP:0009778','MP:0000265','MP:0002883') and phenotype_b.phenotype_id in ('ZP:0006005','ZP:0100399','ZP:0005569','ZP:0137444','ZP:0142224') "
           # "order by phenotype_a.phenotype_id, phenotype_b.phenotype_id").show(max_width=10000, max_rows=50)
           # "where phenotype_a.in_taxon = 'NCBITaxon:10090' and phenotype_b.in_taxon = 'NCBITaxon:9606' ") # .show(max_width=10000, max_rows=50)

print('Phenotype-Ortholog Count Table: ')
duckdb.sql("SELECT * FROM phenotype_ortholog_counts").show(max_width=10000, max_rows=10)

print("Phenotype to ortholog with 0 orthologs (should be empty):")
duckdb.sql("SELECT * FROM phenotype_ortholog_counts where ortholog_count = 0 ").show(max_width=10000, max_rows=10)

# save to disk
duckdb.sql("COPY phenotype_ortholog_counts TO '../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_counts.csv' (HEADER true, DELIMITER '\t');")





# generate common ortholog counts
'''
Here we need to put together a count of distinct ortholog matches per phenotype pair
Can utilize the phenotype_to_ortholog table.
Get the cross-product of all phenotypes, include every ortholog combination between the two phenotypes, have an indicator column for when the orhthologs match
Have a follow up table that counts the 1s. Maybe start with the phenotype table to allow for limiting the inputs?
phenotype_ortholog Columns: phenotype_id	phenotype_name	phenotype_taxon_id	ortholog_id
phenotype columns: phenotype_id	phenotype_name	in_taxon	phenotype_row_number
'''
# Trim phenotype table to remove phenotypes without orthologs:
duckdb.sql("CREATE TABLE tuncated_phenotypes AS "
           "SELECT distinct p.* from phenotypes p "
           "inner join phenotype_ortholog_counts po "
           "on p.phenotype_id = po.phenotype_id "
           "and p.in_taxon = po.phenotype_taxon_id ")

# ERROR: There's a bug here where we're assembling common orthologs between phenotypes without consideration of taxon.
# An example is for HP:0005107 and XPO:0115436 , which should only have rows for XPO:0115436 is in taxon 8364, but we're getting rows for taxon 8355 as well.

# This has high memory usage, longer compute time. Any ways to optimize?
duckdb.sql("CREATE TABLE common_orthologs AS "
           "SELECT distinct p_a.phenotype_id as phenotype_a_id, p_a.in_taxon as phenotype_a_taxon_id, p_b.phenotype_id as phenotype_b_id, p_b.in_taxon as phenotype_b_taxon_id, "
           "po_a.ortholog_id as ortholog_id "
           "FROM tuncated_phenotypes as p_a "
           "LEFT JOIN tuncated_phenotypes as p_b "
           "ON p_a.in_taxon <> p_b.in_taxon " # This will give a cross-product of phenotypes where taxons are not the same.
           "LEFT JOIN phenotype_to_ortholog as po_a "
           "ON p_a.phenotype_id = po_a.phenotype_id and p_a.in_taxon = po_a.phenotype_taxon_id "
           "LEFT JOIN phenotype_to_ortholog as po_b "
           "ON p_b.phenotype_id = po_b.phenotype_id and p_b.in_taxon = po_b.phenotype_taxon_id "
           "WHERE po_a.ortholog_id = po_b.ortholog_id")

print('Common Orthologs Table')
duckdb.sql("SELECT * FROM common_orthologs ").show(max_width=10000, max_rows=10)

print("Common ortholgs: XPO:0115436 should not have any rows for taxon 8355")
duckdb.sql("SELECT * FROM common_orthologs "
           "WHERE phenotype_a_id in ('HP:0005107', 'XPO:0115436' ) "
           "and phenotype_b_id in ('HP:0005107', 'XPO:0115436') "
           "order by phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id").show(max_width=10000, max_rows=100)


# save to disk
duckdb.sql("COPY common_orthologs TO '../../../datasets/intermediate/duckdb_tables/common_orthologs.csv' (HEADER true, DELIMITER '\t');")



duckdb.sql("CREATE TABLE common_ortholog_counts AS "
           "SELECT distinct phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id, common_ortholog_count "
           "FROM (SELECT phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id, count(distinct ortholog_id) as common_ortholog_count "
           "FROM common_orthologs "
           "GROUP BY phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id "
           "UNION ALL "
           "SELECT phenotype_b_id as phenotype_a_id, phenotype_b_taxon_id as phenotype_a_taxon_id, phenotype_a_id as phenotype_b_id, phenotype_a_taxon_id as phenotype_b_taxon_id, count(distinct ortholog_id) as common_ortholog_count "
           "FROM common_orthologs "
           "GROUP BY phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id) ")
           # "where phenotype_a.phenotype_id in ('MP:0001399','MP:0001025','MP:0009778','MP:0000265','MP:0002883') and phenotype_b.phenotype_id in ('ZP:0006005','ZP:0100399','ZP:0005569','ZP:0137444','ZP:0142224') "
           # "order by phenotype_a.phenotype_id, phenotype_b.phenotype_id").show(max_width=10000, max_rows=50)
           # "where phenotype_a.in_taxon = 'NCBITaxon:10090' and phenotype_b.in_taxon = 'NCBITaxon:9606' ") # .show(max_width=10000, max_rows=50)

print('Common Ortholog Count Table')
duckdb.sql("SELECT * FROM common_ortholog_counts ").show(max_width=10000, max_rows=10)

# save to disk
duckdb.sql("COPY common_ortholog_counts TO '../../../datasets/intermediate/duckdb_tables/common_ortholog_counts.csv' (HEADER true, DELIMITER '\t');")




'''
duckdb.sql("CREATE TABLE common_ortholog_counts AS "
           "SELECT  po_a.phenotype_id as phenotype_a_id, po_a.phenotype_taxon_id as phenotype_a_taxon_id, po_b.phenotype_id as phenotype_b_id, po_b.phenotype_taxon_id as phenotype_b_taxon_id, count(distinct po_b.ortholog_id) as common_ortholog_count "
           "FROM phenotype_ortholog as po_a "
           "INNER JOIN phenotype_ortholog as po_b "
           "ON po_a.ortholog_id = po_b.ortholog_id and po_a.phenotype_taxon_id <> po_b.phenotype_taxon_id "
           "GROUP BY po_a.phenotype_id, po_a.phenotype_taxon_id, po_b.phenotype_id, po_b.phenotype_taxon_id "
           "UNION ALL "
           "SELECT  po_b.phenotype_id as phenotype_a_id, po_b.phenotype_taxon_id as phenotype_a_taxon_id, po_a.phenotype_id as phenotype_b_id, po_a.phenotype_taxon_id as phenotype_b_taxon_id, count(distinct po_b.ortholog_id) as common_ortholog_count "
           "FROM phenotype_ortholog as po_a "
           "INNER JOIN phenotype_ortholog as po_b "
           "ON po_a.ortholog_id = po_b.ortholog_id and po_a.phenotype_taxon_id <> po_b.phenotype_taxon_id "
           "GROUP BY po_a.phenotype_id, po_a.phenotype_taxon_id, po_b.phenotype_id, po_b.phenotype_taxon_id ")
           # "where phenotype_a.phenotype_id in ('MP:0001399','MP:0001025','MP:0009778','MP:0000265','MP:0002883') and phenotype_b.phenotype_id in ('ZP:0006005','ZP:0100399','ZP:0005569','ZP:0137444','ZP:0142224') "
           # "order by phenotype_a.phenotype_id, phenotype_b.phenotype_id").show(max_width=10000, max_rows=50)
           # "where phenotype_a.in_taxon = 'NCBITaxon:10090' and phenotype_b.in_taxon = 'NCBITaxon:9606' ") # .show(max_width=10000, max_rows=50)
'''


print('Rows with common_ortholog_count = 0')
duckdb.sql("SELECT * FROM common_ortholog_counts "
           "where common_ortholog_count = 0 ").show(max_width=10000, max_rows=10)



'''
duckdb.sql("CREATE TABLE orthologs AS "
           "select distinct gene_a_id, gene_a_taxon_id, predicate, gene_b_id, gene_b_taxon_id, ortholog_id from ( "
           "select distinct gene_a.gene_id as gene_a_id, gene_A.gene_taxon_id as gene_a_taxon_id, predicate, gene_b.gene_id as gene_b_id, gene_b.gene_taxon_id as gene_b_taxon_id, regexp_replace(has_evidence, '.+:', '')  as ortholog_id "
           "from edges inner join genes as gene_A on edges.subject = gene_a.gene_id "
           "inner join genes as gene_B on edges.object = gene_B.gene_id "
           "where edges.predicate = 'biolink:orthologous_to' "
           "union all "
           "select distinct gene_b.gene_id as gene_a_id, gene_b.gene_taxon_id as gene_a_taxon_id, predicate, gene_a.gene_id as gene_b_id, gene_a.gene_taxon_id as gene_b_taxon_id, regexp_replace(has_evidence, '.+:', '')  as ortholog_id "
           "from edges inner join genes as gene_A on edges.subject = gene_a.gene_id "
           "inner join genes as gene_B on edges.object = gene_B.gene_id "
           "where edges.predicate = 'biolink:orthologous_to' )")
'''






# duckdb.sql("SELECT * FROM phenotype_ortholog_counts").show(max_width=10000, max_rows=10)


# datasets/intermediate/duckdb_tables
# "../../datasets/intermediate/panther/panther_orthologs.tsv"
# datasets/sources/monarch_kg/monarch-kg_edges.tsv





