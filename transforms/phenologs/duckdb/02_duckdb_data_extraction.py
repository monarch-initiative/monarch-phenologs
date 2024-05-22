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
           "and nodes.category = 'biolink:Gene' ")

duckdb.sql("SELECT count(*) as gene_row_count FROM genes").show(max_width=10000, max_rows=100)

duckdb.sql("COPY genes TO '../../../datasets/intermediate/duckdb_tables/genes.csv' (HEADER true, DELIMITER '\t');")



# create phenotypes table
# NOTE: Should this be amended to include the taxon, so that we definitely say phenotype X is in species Y for the rat/mouse phenotypes?
# NOTE: This approach will only include phenotypes that are a phenotype of a gene, so it will drop out phenotypes with no known causal genes.
duckdb.sql("CREATE TABLE phenotypes AS "
           "select distinct nodes_B.id as phenotype_id, "
           "nodes_B.name as phenotype_name, "
           "nodes_A.in_taxon " # This will allow us to assign a taxon for the ambiguous mouse/rat phenotypes   
           "from nodes as nodes_A left join edges on nodes_A.id = edges.subject "
           "left join nodes as nodes_B "
           "on edges.object = nodes_B.id "
           "where nodes_A.category = 'biolink:Gene' "
           "and edges.predicate = 'biolink:has_phenotype' "
           "and nodes_B.category = 'biolink:PhenotypicFeature' "     
           "and (regexp_replace(nodes_B.id, ':.+', ':') like 'HP:%' "
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'MP:%' " # Both rat and mouse use MP
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'ZP:%' "
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'XPO:%' "
           "or regexp_replace(nodes_B.id, ':.+', ':') like 'WBPhenotype:%' "
           ")")



duckdb.sql("SELECT count(*) as phenotype_row_count FROM phenotypes").show(max_width=10000, max_rows=100)

duckdb.sql("COPY phenotypes TO '../../../datasets/intermediate/duckdb_tables/phenotypes.csv' (HEADER true, DELIMITER '\t');")




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


duckdb.sql("SELECT * FROM orthologs "
           "where (gene_a_id = 'HGNC:16435' and gene_b_id = 'MGI:95886') or "
           "(gene_b_id = 'HGNC:16435' and gene_a_id = 'MGI:95886') ").show(max_width=10000, max_rows=10)

duckdb.sql("SELECT count(*) as ortholog_row_count FROM orthologs").show(max_width=10000, max_rows=10)
# total count: 178051
# distinct count: 178030
# bidirectional via union all: 356060

# write to file
duckdb.sql("COPY orthologs TO '../../../datasets/intermediate/duckdb_tables/orthologs.csv' (HEADER true, DELIMITER '\t');")




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

# write to file
duckdb.sql("COPY gene_to_phenotype TO '../../../datasets/intermediate/duckdb_tables/gene_to_phenotype.csv' (HEADER true, DELIMITER '\t');")




# create gene to ortholog table, subset to common orthologs table for taxons of interest?


# create phenotype-ortholog file
# gene to phenotype schema: gene_id	gene_name	gene_taxon_id	gene_taxon_label	predicate	phenotype_id	phenotype_name
# ortholog schema: gene_a_id	gene_a_taxon_id	predicate	gene_b_id	gene_b_taxon_id	ortholog_id
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS "
           "SELECT distinct g2p.phenotype_id, g2p.phenotype_name, g2p.phenotype_taxon_id, orth.ortholog_id FROM gene_to_phenotype as g2p left join orthologs as orth on g2p.gene_id = orth.gene_a_id ")

# write to file
duckdb.sql("COPY phenotype_to_ortholog TO '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv' (HEADER true, DELIMITER '\t');")

duckdb.sql("SELECT count(*) as phenotype_ortholog_row_count FROM phenotype_to_ortholog").show(max_width=10000, max_rows=10)
# total count: 1999745
# distinct count: 447587 -> after adding phenotype taxon ID: 448251 -> Without the phenotype taxon, the overlapping mouse/rat phenotypes collapse.




# datasets/intermediate/duckdb_tables
# "../../datasets/intermediate/panther/panther_orthologs.tsv"
# datasets/sources/monarch_kg/monarch-kg_edges.tsv





