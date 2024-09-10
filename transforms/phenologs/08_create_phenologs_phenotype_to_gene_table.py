import pandas as pd
import pickle
import csv
import duckdb




'''
As a final step, generate a table from the significant phenologs datafile that includes a list of genes on each side of the phenolog relationship.

Could have a list of genes for each species, list of list of common orthologs, perhaps a separate list that indicates which genes were the common genes/orthologs between the two species?

Columns for current significant phenologs file: 
Species A Phenotype ID	Species B Phenotype ID	p-value	Significance


Output table format:
Species A Phenotype ID, Species A Gene ID List, Species A Gene Name List, Species B Phenotype ID, Species B Gene ID List, Species B Gene Name List, p-value, Significance
'''
duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Create the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../datasets/intermediate/duckdb_tables/genes.csv'") # Create the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Create the orthologs table
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS SELECT * FROM '../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv'") # Create the phenotype to ortholog table
duckdb.sql("CREATE TABLE phenolog_base AS SELECT * FROM '../../datasets/intermediate/duckdb_tables/phenolog_base.csv'") # Create the phenolog base table
duckdb.sql("CREATE TABLE sig_phenologs AS SELECT * FROM '../../datasets/output/phenologs/significant_phenolog_data.tsv'") # Load the significant phenologs table

duckdb.sql("SELECT * FROM sig_phenologs").show(max_width=10000, max_rows=10)

# TODO: This process will need to be udpated to account for the phenotaxon IDs and not just the phenotypes.
# Should preferably just update the output file generation.
# Create a phenolog to orhtolog table that assembles arrays of orthologs for each phenotpe in the phenolog pair.

duckdb.sql("ALTER TABLE sig_phenologs RENAME 'Species A Phenotype ID' TO 'species_a_phenotype_id'")
duckdb.sql("ALTER TABLE sig_phenologs RENAME 'Species B Phenotype ID' TO 'species_b_phenotype_id'")
duckdb.sql("ALTER TABLE sig_phenologs RENAME 'p-value' TO 'p_value'")
duckdb.sql("ALTER TABLE sig_phenologs RENAME 'Significance' TO 'significance'")

duckdb.sql("SELECT * FROM sig_phenologs").show(max_width=10000, max_rows=10)

# Helper table: A pre-calculated table of ortholog overlaps between phenotypes.
# This could probably be built from the phenolog_base table
duckdb.sql("CREATE TABLE phenotype_ortholog_overlap AS "
           "SELECT * FROM phenolog_base pa "
           "LEFT JOIN phenotypes pb "
           "on "
           "WHERE phenotype_a = 'HP:0008247' and phenotype_b = 'MP:0005586 ' ")





duckdb.sql("SELECT * FROM phenotype_ortholog_overlap").show(max_width=10000, max_rows=10)


'''
duckdb.sql("CREATE TABLE phenolog_to_ortholog AS "
           "SELECT "
           "s.species_a_phenotype_id, s.species_b_phenotype_id, s.p_value, s.significance, "
           ""
           "FROM sig_phenologs s "
           "LEFT JOIN phenotypes pa "
           "ON s.species_a_phenotype_id = pa.phenotype_id "
           "LEFT JOIN phenotypes pb "
           "ON s.species_b_phenotype_id = pb.phenotype_id "
           "LEFT JOIN phenotype_to_ortholog poa "
           "ON s.species_a_phenotype_id = poa.phenotype_id "
           "LEFT JOIN phenotype_to_ortholog pob "
           "ON s.species_b_phenotype_id = pob.phenotype_id "
           "WHERE s.species_a_phenotype_id = 'HP:0008247' "
           "GROUP BY s.species_a_phenotype_id, s.species_b_phenotype_id, s.p_value, s.significance ")


duckdb.sql("SELECT * FROM phenolog_to_ortholog").show(max_width=10000, max_rows=10)

'''