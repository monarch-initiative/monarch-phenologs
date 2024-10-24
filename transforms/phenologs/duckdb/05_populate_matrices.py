import duckdb
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, pearsonr
from duckdb.typing import *


'''
If we already have the phenolog base table with hypergeometric probability/weight matrix, 
let's fill in the distance matrix using Pearson correlation. 

'''




# Load tables
duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Create the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/genes.csv'") # Create the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Create the genes table
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv'") # Create the phenotype-ortholog table

duckdb.sql("CREATE TABLE common_orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_orthologs.csv'") # Load the common orthologs table
duckdb.sql("CREATE TABLE common_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_ortholog_counts.csv'") # Load the common ortholog counts table
duckdb.sql("CREATE TABLE phenotype_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_counts.csv'") # Load the phenotype-ortholog counts table

duckdb.sql("CREATE TABLE phenolog_base AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/updated_phenolog_base.csv'") # Load the phenologs base table


# Need to generate a list of taxons for reference:
duckdb.sql("CREATE TABLE taxons as SELECT distinct phenotype_taxon_id as taxon_id FROM phenotype_to_ortholog")
print('Taxons Table: ')
duckdb.sql("SELECT * FROM taxons").show(max_width=10000, max_rows=10)
# Use the common orthologs table for the ortholog portion of the matrix or the base orthologs table filtered to included taxons?
# TODO: I'm thinking the latter is more appropriate, but double-check.
# Also, should we be using row numbers - 1 to align with the matrix coordinates?
duckdb.sql("CREATE TABLE knn_orthologs AS "
           "SELECT row_number() OVER (ORDER BY (ortholog_id)) - 1 as j, * FROM ( "
           "SELECT DISTINCT ortholog_id from orthologs "
           "WHERE "
           "gene_a_taxon_id in (select taxon_id from taxons) "
           "and "
           "gene_b_taxon_id in (select taxon_id from taxons) "
           "order by ortholog_id )")

print('kNN Orthologs Table: ')
duckdb.sql("SELECT * FROM knn_orthologs").show(max_width=10000, max_rows=10)

print('kNN Orthologs Table row count: ')
duckdb.sql("SELECT count(*) FROM knn_orthologs").show(max_width=10000, max_rows=10)

# Need to generate a phenotype +10 nearest neighbors table
# *****Create phenotype IDs as phenotype_taxon?
# Might instead need to create an array column instead of 10 distinct columns, allowing for the situation where 10+ nearest neighbors are found when distance is the same.
# TODO: However, this could result in there be many more than 10 nearest neighbors. Should run a data check and see if there's any major outliers.
duckdb.sql("CREATE TABLE knn_phenotypes AS "
           "SELECT row_number() OVER (ORDER BY (phenotype_id))- 1 as i, * FROM ( "
           "SELECT DISTINCT pto.phenotype_id, pto.phenotype_name, pto.phenotype_taxon_id, concat(pto.phenotype_id, '_', pto.phenotype_taxon_id) as phenotaxon_id, "
           # "row_number() OVER (ORDER BY (SELECT NULL)) as j, "
           "null as knn_phenotype_array, "
           "null as knn_phenotype_1, "
           "null as knn_taxon_1, "
           "null as knn_phenotype_2, "
           "null as knn_taxon_2, "
           "null as knn_phenotype_3, "
           "null as knn_taxon_3, "
           "null as knn_phenotype_4, "
           "null as knn_taxon_4, "
           "null as knn_phenotype_5, "
           "null as knn_taxon_5, "
           "null as knn_phenotype_6, "
           "null as knn_taxon_6, "
           "null as knn_phenotype_7, "
           "null as knn_taxon_7, "
           "null as knn_phenotype_8, "
           "null as knn_taxon_8, "
           "null as knn_phenotype_9, "
           "null as knn_taxon_9, "
           "null as knn_phenotype_10, "
           "null as knn_taxon_10 "
           "FROM phenotype_to_ortholog pto )")
           # "order by j")


print('kNN Phenotypes Table: ')
duckdb.sql("SELECT * FROM knn_phenotypes order by phenotaxon_id").show(max_width=10000, max_rows=10)

print('kNN Phenotypes Table row count: ')
duckdb.sql("SELECT count(*) FROM knn_phenotypes").show(max_width=10000, max_rows=10)




# In preparation for filling out the matrix, generate the phenotype_to_ortholog table with the phenotaxon_id added:
duckdb.sql("CREATE TABLE phenotype_ortholog_matrix_input AS "
           "SELECT DISTINCT pto.phenotype_id, pto.phenotype_name, pto.phenotype_taxon_id, "
           "concat(pto.phenotype_id, '_', pto.phenotype_taxon_id) as phenotaxon_id, pto.ortholog_id, "
           "p.i, o.j "
           "FROM phenotype_to_ortholog pto "
           "LEFT JOIN knn_phenotypes p "
           "ON concat(pto.phenotype_id, '_', pto.phenotype_taxon_id) = p.phenotaxon_id "
           "LEFT JOIN knn_orthologs o "
           "ON pto.ortholog_id = o.ortholog_id")

print('phenotype_ortholog_matrix_input Table: ')
duckdb.sql("SELECT * FROM phenotype_ortholog_matrix_input order by phenotaxon_id").show(max_width=10000, max_rows=10)

print('phenotype_ortholog_matrix_input row count: ')
duckdb.sql("SELECT count(*) FROM phenotype_ortholog_matrix_input").show(max_width=10000, max_rows=10)



# Need to generate a phenotype-ortholog matrix that can be referred to/sliced in later operations
# duckdb.sql("CREATE TABLE knn_phenotypes AS SELECT DISTINCT phenotype_id, phenotype_taxon_id from phenotype_to_ortholog")
# utilize the knn_phenotypes table, which has the merged phenotaxon column.
# Create a zero matrix of length knn_phenotypes x knn_orthologs

ortholog_count = duckdb.sql("SELECT count(*) FROM knn_orthologs").fetchall()
ortholog_count = ortholog_count[0][0]
print("Ortholog count:" + str(ortholog_count))

phenotype_count = duckdb.sql("SELECT count(*) FROM knn_phenotypes").fetchall()
phenotype_count = phenotype_count[0][0]
print("Phenotype count:" + str(phenotype_count))

matrix_size = phenotype_count * ortholog_count
print("Phenotype-ortholog matrix size (total cells): " + str(matrix_size))

# pandas test
# pandas_count = duckdb.df("SELECT count(*) FROM knn_phenotypes").fetchall()




# So just create a matrix of size # phenotaxons x # orthologs
# If desired, can number the phenotaxons and orhtologs to use as indexes: phenotype i, ortholog j

# We have the phenotypes and the orthologs available with numeric identifiers/row numbers starting from 0,
# and these could be used as phenotype-ortholog matrix coordinates.

phenotype_ortholog_matrix = np.zeros((phenotype_count, ortholog_count), dtype=np.int8)

'''
FOR REFERENCE: Phenotypes: i, Orthologs: j
Conceptualizing the phenotype ortholog matrix and the downstream calculations:
Note that the phenotype-ortholog matrix is a matrix of all orthologs that are present for the included taxons, 
regardless of whether or not an association of that ortholog with a phenotype exists. 
In other words, we're including the full "ortholog space."
Then when we compare two phenotypes, we take the representative slices of the matrix for each phenotype
(or more specifically each "phenotaxon") and compare those slices using Pearson correlation.
'''



# Here need to iterate through and add a 1 for every row in knn_phenotypes where there's a corresponding row in knn
# What's the best approach? Convert the duckdb tables to pandas/polars dataframes, iterate through?

# STEPS:
# Iterate through the rows in the newly constructed phenotype_ortholog_matrix_input table:
# Columns: phenotype_id, phenotype_name, phenotype_taxon_id, phenotaxon_id, ortholog_id
# Self comparison? Where phenotaxon_id <> phenotaxon_id and phenotype_taxon_id <> phenotype_taxon_id?

# Place a 1 at the matrix coordinates corresponding to the ith phenotype, jth ortholog.
# Best way of extracting these coordinate values? Should add the i/j columns to the phenotype_ortholog_matrix_input


# read the result of an arbitrary SQL query to a Pandas DataFrame
coordinates_df = duckdb.sql("SELECT i, j from phenotype_ortholog_matrix_input").df()
print("coordinates_df:")
print(coordinates_df)
# counter = 0
for row in coordinates_df.itertuples():
    # print(row)
    i = row.i
    j = row.j
    # counter += 1
    # print("i: " + str(i) + ", j: " + str(j))
    phenotype_ortholog_matrix[i][j] = 1

print(phenotype_ortholog_matrix[10333][3368])

np.save('../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_matrix.npy', phenotype_ortholog_matrix)
np.savetxt('../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_matrix.txt', phenotype_ortholog_matrix)

# Just to test, take two rows of the phenotype_ortholog_matrix and compare:
row_a = phenotype_ortholog_matrix[10333]
row_b = phenotype_ortholog_matrix[10332]
print(row_a)
print(row_b)
(coefficient, p_value) = pearsonr(phenotype_ortholog_matrix[10333],
                                          phenotype_ortholog_matrix[10332])

print("Forward: ")
print(coefficient)
print(p_value)

(coefficient, p_value) = pearsonr(phenotype_ortholog_matrix[10333],
                                          phenotype_ortholog_matrix[10332])
print("Reverse:")
print(coefficient)
print(p_value)
# TODO:
# Would need to iterate through every phenotaxon_id vs phenotaxon_id (aka updated phenologs base table)
# and grab the pearson coefficient for the corresponding rows of the phenotype-ortholog matrix.

# NOTE: May need to move this to a separate script to avoid overburdening the runtime.

# What would work best? Would it work to pass data out of DuckDB, or would it work better to create a functionw within DuckDB as done elsewhere?




''' DuckDB function to calculate pearson correlation:
In order to calculate the pearson correlation we need to pass along the phenotaxon IDs of both phenotypes.
We also need to refer the phenotype_ortholog_matrix and extract row slices -> Any issue referencing that matrix from within a DuckDB query? Let's find out....

'''

def calculate_pearson_correlation_coefficient(phenotaxon_a: int, phenotaxon_b: int):
    # Need to convert the phenotaxon IDs to i/j coordinates? These should be available within the phenologs base table


    (coefficient, p_value) = pearsonr(phenotype_ortholog_matrix[phenotaxon_a],
                                      phenotype_ortholog_matrix[phenotaxon_b])


    # print("Phenotype B Ortholog Count: " + str(phenotype_b_ortholog_count) + ", Phenotype A Ortholog Count: " + str(phenotype_a_ortholog_count) + ", Shared ortholog count: "+ str(shared_ortholog_count) + ", Ortholog matches: " + str(ortholog_matches) + ", Probability: " + str(probability))
    # return (coefficient, p_value)
    return coefficient

def calculate_pearson_correlation_p_value(phenotaxon_a: int, phenotaxon_b: int):
    # Need to convert the phenotaxon IDs to i/j coordinates? These should be available within the phenologs base table


    (coefficient, p_value) = pearsonr(phenotype_ortholog_matrix[phenotaxon_a],
                                      phenotype_ortholog_matrix[phenotaxon_b])


    # print("Phenotype B Ortholog Count: " + str(phenotype_b_ortholog_count) + ", Phenotype A Ortholog Count: " + str(phenotype_a_ortholog_count) + ", Shared ortholog count: "+ str(shared_ortholog_count) + ", Ortholog matches: " + str(ortholog_matches) + ", Probability: " + str(probability))
    # return (coefficient, p_value)
    return p_value

duckdb.create_function("pearson_corr_coeff", calculate_pearson_correlation_coefficient, [ int, int], DOUBLE)

duckdb.create_function("pearson_corr_p_value", calculate_pearson_correlation_p_value, [ int, int], DOUBLE)


test = calculate_pearson_correlation_coefficient(10333, 10332)

print("Function test - coefficient: ")
print(test)


test = calculate_pearson_correlation_p_value(10333, 10332)

print("Function test - p-value: ")
print(test)

print("phenologs base table: ")
duckdb.sql("SELECT * FROM phenolog_base ").show(max_width=10000, max_rows=10)

# Limit table size for testing.
# duckdb.sql("REPLACE TABLE phenolog_base AS SELECT * from phenolog_base limit 1000")

duckdb.sql("ALTER TABLE phenolog_base ADD COLUMN pearson_coefficient DOUBLE;")
duckdb.sql("ALTER TABLE phenolog_base ADD COLUMN pearson_p_value DOUBLE;")
duckdb.sql("ALTER TABLE phenolog_base ADD COLUMN phenotaxon_a_matrix_coordinate INTEGER;")
duckdb.sql("ALTER TABLE phenolog_base ADD COLUMN phenotaxon_b_matrix_coordinate INTEGER;")

print("Updated phenologs base table: ")
duckdb.sql("SELECT * FROM phenolog_base ").show(max_width=10000, max_rows=10)






# Set the phenotaxon coordinates from the
print("Updating phenotaxon A coordinates.")
duckdb.sql("UPDATE phenolog_base "
           "SET phenotaxon_a_matrix_coordinate = x.i "
           "FROM knn_phenotypes x "
           "WHERE phenotaxon_a_id = x.phenotaxon_id ")

print("Updating phenotaxon B coordinates.")
duckdb.sql("UPDATE phenolog_base "
           "SET phenotaxon_b_matrix_coordinate = x.i "
           "FROM knn_phenotypes x "
           "WHERE phenotaxon_b_id = x.phenotaxon_id ")

print("Updated phenologs base table: ")
duckdb.sql("SELECT * FROM phenolog_base ").show(max_width=10000, max_rows=10)


# Possible to update two columns at once so we don't have to run through the entire table twice?
print("Performing pearson correlation calculations.")
duckdb.sql("UPDATE phenolog_base "
           "SET pearson_coefficient = pearson_corr_coeff(phenotaxon_a_matrix_coordinate, phenotaxon_b_matrix_coordinate), pearson_p_value = pearson_corr_p_value(phenotaxon_a_matrix_coordinate, phenotaxon_b_matrix_coordinate) ")


print("Updated phenologs base table: populated p-values?")
duckdb.sql("SELECT * FROM phenolog_base ").show(max_width=10000, max_rows=10)


# save to disk -> Overwrite?
duckdb.sql("COPY phenolog_base TO '../../../datasets/intermediate/duckdb_tables/05_updated_phenolog_base.csv' (HEADER true, DELIMITER '\t');")


duckdb.sql("COPY phenotype_ortholog_matrix_input TO '../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_matrix_input.csv' (HEADER true, DELIMITER '\t');")
duckdb.sql("COPY knn_orthologs TO '../../../datasets/intermediate/duckdb_tables/knn_orthologs.csv' (HEADER true, DELIMITER '\t');")
duckdb.sql("COPY knn_phenotypes TO '../../../datasets/intermediate/duckdb_tables/knn_phenotypes.csv' (HEADER true, DELIMITER '\t');")
duckdb.sql("COPY taxons TO '../../../datasets/intermediate/duckdb_tables/taxons.csv' (HEADER true, DELIMITER '\t');")
