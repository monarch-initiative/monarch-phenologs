import duckdb
import os

duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Create the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/genes.csv'") # Create the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Create the genes table
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv'") # Create the phenotype-ortholog table

duckdb.sql("CREATE TABLE common_orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_orthologs.csv'") # Load the common orthologs table
duckdb.sql("CREATE TABLE common_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_ortholog_counts.csv'") # Load the common ortholog counts table
duckdb.sql("CREATE TABLE phenotype_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_counts.csv'") # Load the phenotype-ortholog counts table







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
# Necessary to stitch together the various taxon pair-wise combinations?'
# Can we trim down the inputs? How about we drop all phenotypes not present in the phenotype-ortholog table?

# Trim phenotype table to remove phenotypes without orthologs:
duckdb.sql("CREATE TABLE tuncated_phenotypes AS "
           "SELECT distinct p.* from phenotypes p "
           "inner join phenotype_ortholog_counts po "
           "on p.phenotype_id = po.phenotype_id "
           "and p.in_taxon = po.phenotype_taxon_id ")
# phenotype_id	phenotype_name	in_taxon	phenotype_row_number
# phenotype_id	phenotype_name	phenotype_taxon_id	ortholog_id

print("Phenotype row count: ")
duckdb.sql("SELECT count(*) as row_count FROM phenotypes").show(max_width=10000, max_rows=10)
print("Truncated phenotype row count: ")
duckdb.sql("SELECT count(*) as row_count FROM tuncated_phenotypes").show(max_width=10000, max_rows=10)


# A few different performance tuning options:
# duckdb.sql("SET threads = 20") # This may work as setting a limit of threads which may be less than what DuckDB would use normally.
# duckdb.sql("PRAGMA max_temp_directory_size='100GiB'")
duckdb.sql("CREATE TABLE phenolog_base AS "
           "SELECT DISTINCT pa.phenotype_id as phenotype_a_id, pa.in_taxon as phenotype_a_in_taxon, pb.phenotype_id as phenotype_b_id, pb.in_taxon as phenotype_b_in_taxon "
           # ", null as phenotype_a_ortholog_count, null as phenotype_b_ortholog_count, null as shared_ortholog_count, null as ortholog_matches, "
           # "null as hypergeom_probability, "
           # "null as p_value, null as significance, "
           # "null as distance, null as weight "
           "FROM tuncated_phenotypes pa "
           "LEFT JOIN tuncated_phenotypes pb "
           "ON pa.phenotype_id <> pb.phenotype_id and pa.in_taxon <> pb.in_taxon ")
           # "where phenotype_a.phenotype_id in ('MP:0001399','MP:0001025','MP:0009778','MP:0000265','MP:0002883') and phenotype_b.phenotype_id in ('ZP:0006005','ZP:0100399','ZP:0005569','ZP:0137444','ZP:0142224') "
           # "order by phenotype_a.phenotype_id, phenotype_b.phenotype_id").show(max_width=10000, max_rows=50)
           # "where phenotype_a.in_taxon = 'NCBITaxon:10090' and phenotype_b.in_taxon = 'NCBITaxon:9606' ") # .show(max_width=10000, max_rows=50)

print('Phenolog base table: ')
duckdb.sql("SELECT * FROM phenolog_base").show(max_width=10000, max_rows=10)


print('Phenolog base table row count: ')
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