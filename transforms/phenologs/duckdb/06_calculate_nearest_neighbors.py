import duckdb
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, pearsonr
from duckdb.typing import *


'''
Purpose:
Need to identify the top 10 phenotaxons for each phenotaxon where the species are not the same.
Could do this in SQL by using a row number by partition function.
'''


# Load tables
duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Create the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/genes.csv'") # Create the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Create the genes table
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv'") # Create the phenotype-ortholog table

duckdb.sql("CREATE TABLE common_orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_orthologs.csv'") # Load the common orthologs table
duckdb.sql("CREATE TABLE common_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_ortholog_counts.csv'") # Load the common ortholog counts table
duckdb.sql("CREATE TABLE phenotype_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_counts.csv'") # Load the phenotype-ortholog counts table

#TODO: Update to 05_ table version when available.
duckdb.sql("CREATE TABLE phenolog_base AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/updated_phenolog_base.csv'") # Load the phenologs base table

duckdb.sql("SELECT * FROM phenolog_base "
           "WHERE (phenotaxon_a_id = 'MP:0006277_NCBITaxon:10090' and phenotaxon_b_id = 'HP:0100856_NCBITaxon:9606') "
           "or (phenotaxon_a_id = 'HP:0100856_NCBITaxon:9606' and phenotaxon_b_id = 'MP:0006277_NCBITaxon:10090')").show(max_width=10000, max_rows=10)

# Could try doing a row number over partition and see if that could quickly identify the k-nearest neighbor rows.
# Available columns:
# phenotype_a_id	phenotype_a_in_taxon	phenotaxon_a_id	phenotype_b_id	phenotype_b_in_taxon	phenotaxon_b_id	phenotype_a_ortholog_count	phenotype_b_ortholog_count
# ortholog_matches	shared_orthologs	p_value	pearson_coefficient	pearson_p_value	phenotaxon_a_matrix_coordinate	phenotaxon_b_matrix_coordinate
# May have to get creative since instances of all phenotype pairs for a given phenotype could have that phenotype in column a or column b.
'''
Possible approach: 
-Union the table with itself and the phenotype slots flipped.
-Do the row over partition when grouped by phenotaxon A.
-Grab the top 10 rows for phenotaxon A
-Add those phenotaxon IDs to an array column in the phenotypes (?) table?
'''



duckdb.sql("CREATE TABLE knn_table AS "
           "SELECT phenotype_a_id, phenotype_a_in_taxon, phenotaxon_a_id, phenotype_b_id, phenotype_b_in_taxon, phenotaxon_b_id, pearson_coefficient, pearson_p_value "
           "FROM phenolog_base "
           "UNION ALL "
           "SELECT phenotype_b_id as phenotype_a_id, phenotype_b_in_taxon as phenotype_a_in_taxon, phenotaxon_b_id as phenotaxon_a_id, "
           "phenotype_a_id as phenotype_b_id, phenotype_a_in_taxon as phenotype_b_in_taxon, phenotaxon_a_id as phenotaxon_b_id, "
           "pearson_coefficient, pearson_p_value"
           "FROM phenolog_base ")

# "row_number() OVER (PARTITION BY list_sort(array_value(concat(pa.phenotype_id, '_', pa.in_taxon), concat(pb.phenotype_id, '_', pb.in_taxon)))) as phenolog_row "
"row_number() OVER (PARTITION BY ORDER BY ) as knn_rank  "

print("kNN Table: ")
duckdb.sql("SELECT * FROM knn_table").show(max_width=10000, max_rows=10)


duckdb.sql("CREATE TABLE knn_rank AS "
           "SELECT *, "
           "row_number() OVER (PARTITION BY phenotaxon_a_id ORDER BY pearson_coefficient) as knn_rank "
           "FROM knn_table ")

print("knn Rank Table: ")
duckdb.sql("SELECT * FROM knn_rank").show(max_width=10000, max_rows=10)

print("Top 10 for phenotaxon ID: MP:0000137_NCBITaxon:10090")
duckdb.sql("SELECT * FROM knn_rank where phenotaxon_a_id = 'MP:0000137_NCBITaxon:10090  '").show(max_width=10000, max_rows=10)
