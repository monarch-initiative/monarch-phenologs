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

duckdb.sql("CREATE TABLE sig_phenologs AS SELECT * FROM '../../datasets/output/phenologs/significant_phenolog_data.tsv'") # Load the significant phenologs table
duckdb.sql("CREATE TABLE sig_phenologs AS SELECT * FROM '../../datasets/output/phenologs/significant_phenolog_data.tsv'") # Load the phenotypes table


duckdb.sql("CREATE TABLE sig_phenologs AS SELECT * FROM '../../datasets/output/phenologs/significant_phenolog_data.tsv'") # Load the genes table
duckdb.sql("SELECT * FROM sig_phenologs").show(max_width=10000, max_rows=10)


