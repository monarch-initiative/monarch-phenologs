import duckdb
import os

duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Create the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/genes.csv'") # Create the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Create the genes table
duckdb.sql("CREATE TABLE phenotype_ortholog AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv'") # Create the phenotype-ortholog table



# generate phenotype-ortholog counts
'''
Here we need to put together a simple count of distinct orthologs per phenotype
Can utilize the phenotype_to_ortholog table.
Columns: phenotype_id	phenotype_name	phenotype_taxon_id	ortholog_id
'''

duckdb.sql("CREATE TABLE phenotype_ortholog_counts AS "
           "SELECT DISTINCT phenotype_id, phenotype_taxon_id, count(distinct ortholog_id) "
           "FROM phenotype_ortholog "
           "GROUP BY phenotype_id, phenotype_taxon_id ")
           # "where phenotype_a.phenotype_id in ('MP:0001399','MP:0001025','MP:0009778','MP:0000265','MP:0002883') and phenotype_b.phenotype_id in ('ZP:0006005','ZP:0100399','ZP:0005569','ZP:0137444','ZP:0142224') "
           # "order by phenotype_a.phenotype_id, phenotype_b.phenotype_id").show(max_width=10000, max_rows=50)
           # "where phenotype_a.in_taxon = 'NCBITaxon:10090' and phenotype_b.in_taxon = 'NCBITaxon:9606' ") # .show(max_width=10000, max_rows=50)

duckdb.sql("SELECT * FROM phenotype_ortholog_counts").show(max_width=10000, max_rows=10)


# generate common ortholog counts
'''
Here we need to put together a count of distinct ortholog matches per phenotype pair
Can utilize the phenotype_to_ortholog table.
Get the cross-product of all phenotypes, include every ortholog combination between the two phenotypes, have an indicator column for when the orhthologs match
Have a follow up table that counts the 1s. Maybe start with the phenotype table to allow for limiting the inputs?
phenotype_ortholog Columns: phenotype_id	phenotype_name	phenotype_taxon_id	ortholog_id
phenotype columns: phenotype_id	phenotype_name	in_taxon	phenotype_row_number
'''

duckdb.sql("CREATE TABLE common_orthologs AS "
           "SELECT distinct p_a.phenotype_id as phenotype_a_id, p_a.in_taxon as phenotype_a_taxon_id, p_b.phenotype_id as phenotype_b_id, p_b.in_taxon as phenotype_b_taxon_id, "
           "po_a.ortholog_id as ortholog_id "
           "FROM phenotypes as p_a "
           "LEFT JOIN phenotypes as p_b "
           "ON p_a.in_taxon <> p_b.in_taxon "
           "LEFT JOIN phenotype_ortholog as po_a "
           "ON p_a.phenotype_id = po_a.phenotype_id "
           "LEFT JOIN phenotype_ortholog as po_b "
           "ON p_b.phenotype_id = po_b.phenotype_id "
           "WHERE po_a.ortholog_id = po_b.ortholog_id")

duckdb.sql("SELECT * FROM common_orthologs ").show(max_width=10000, max_rows=10)

duckdb.sql("CREATE TABLE commong_ortholog_counts AS "
           "SELECT phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id, count(distinct ortholog_id) as common_ortholog_count "
           "FROM common_orthologs "
           "GROUP BY phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id "
           "UNION ALL "
           "SELECT phenotype_b_id as phenotype_a_id, phenotype_b_taxon_id as phenotype_a_taxon_id, phenotype_a_id as phenotype_b_id, phenotype_a_taxon_id as phenotype_b_taxon_id, count(distinct ortholog_id) as common_ortholog_count "
           "FROM common_orthologs "
           "GROUP BY phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id ")
           # "where phenotype_a.phenotype_id in ('MP:0001399','MP:0001025','MP:0009778','MP:0000265','MP:0002883') and phenotype_b.phenotype_id in ('ZP:0006005','ZP:0100399','ZP:0005569','ZP:0137444','ZP:0142224') "
           # "order by phenotype_a.phenotype_id, phenotype_b.phenotype_id").show(max_width=10000, max_rows=50)
           # "where phenotype_a.in_taxon = 'NCBITaxon:10090' and phenotype_b.in_taxon = 'NCBITaxon:9606' ") # .show(max_width=10000, max_rows=50)

duckdb.sql("SELECT * FROM commong_ortholog_counts ").show(max_width=10000, max_rows=10)

'''
duckdb.sql("CREATE TABLE commong_ortholog_counts AS "
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
duckdb.sql("SELECT * FROM commong_ortholog_counts "
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





duckdb.sql("SELECT * FROM phenotype_ortholog_counts").show(max_width=10000, max_rows=10)



# generate the phenolog table

'''
What's needed here?
-A cross-product of every phenotype vs every phenotype where the taxons don't align. 
    NOTE: Since rat uses MP, need to create a subset of MP that actually is found in rat. Same for mouse?
    Worth double-checking: Is there a collision between mouse and rat phenotypes and orthologs? Is it problematic to merge rat and mouse phenotypes into one table?
    Also, don't duplicate the comparisons (e.g. phenotype A vs phenotype B as well as phenotype B vs phenotype A)
-Build onto this table the data needed to perform the hypergeometric probability:
    -Orthologs attached to phenotype A, orthologs attached to phenotype B -> Can basicallly do a row count of the phenotype to orhtholog table grouped by phenotype
    -Orthologs that match between the two phenotypes -> Can basically do a row count on an inner join on orthologs between phenotypes A & B
    -Total common orthologs between the two species -> Can precalculate for all species combos, 

# Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
# N = total number of orthologs shared between species
# n = number of orthologs in species A phenotype -> Can join to phenotype_ortholog_counts table
# m = number of orthologs in species B phenotype -> Can join to phenotype_ortholog_counts table
# c = number of common orthologs between phenotypes (ortholog matches)
'''


# phenotypes schema: phenotype_id	phenotype_name	in_taxon

# Performance here may not be sufficient to create the full cross-product table in one go.
# Necessary to stitch together the various taxon pair-wise combinations?

duckdb.sql("CREATE TABLE phenolog_base AS "
           "SELECT DISTINCT phenotype_a.phenotype_id as phenotype_a_id, phenotype_a.in_taxon, phenotype_b.phenotype_id as phenotype_b_id, phenotype_b.in_taxon, "
           "null as phenotype_a_ortholog_count, null as phenotype_b_ortholog_count, null as shared_ortholog_count, null as ortholog_matches "
           "FROM phenotypes as phenotype_a "
           "LEFT JOIN phenotypes as phenotype_b "
           "ON phenotype_a.phenotype_id <> phenotype_b.phenotype_id and phenotype_a.in_taxon <> phenotype_b.in_taxon ")
           # "where phenotype_a.phenotype_id in ('MP:0001399','MP:0001025','MP:0009778','MP:0000265','MP:0002883') and phenotype_b.phenotype_id in ('ZP:0006005','ZP:0100399','ZP:0005569','ZP:0137444','ZP:0142224') "
           # "order by phenotype_a.phenotype_id, phenotype_b.phenotype_id").show(max_width=10000, max_rows=50)
           # "where phenotype_a.in_taxon = 'NCBITaxon:10090' and phenotype_b.in_taxon = 'NCBITaxon:9606' ") # .show(max_width=10000, max_rows=50)

duckdb.sql("SELECT count(*) as phenolog_base_row_count FROM phenolog_base").show(max_width=10000, max_rows=10)


duckdb.sql("COPY phenolog_base TO '../../../datasets/intermediate/duckdb_tables/phenolog_base.csv' (HEADER true, DELIMITER '\t');")


'''
Next step, want to assemble the variable columns necessary for the hypergeometric probability calculation: 
m = float(phenotype_b_ortholog_count)
n = float(phenotype_a_ortholog_count)
N = float(shared_ortholog_count)
c = float(ortholog_matches)

How to add a column to an existing table? Alter?
'''
# duckdb.sql("ALTER TABLE phenolog_base ADD_SELECT * FROM phenolog_base;")




'''
duckdb.sql("CREATE TABLE phenolog_base AS "
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
'''