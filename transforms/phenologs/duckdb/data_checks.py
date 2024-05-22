# Ideally anything here could be moved into a unit test.

import duckdb

duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Create the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/genes.csv'") # Create the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Create the genes table


# Check that all genes have a taxon
duckdb.sql("select count(*) as empty_gene_taxon_count from genes where gene_taxon_id is null or gene_taxon_id = ''").show()

# Check that all phenotypes have a taxon
duckdb.sql("select count(*) as empty_phenotype_taxon_count from phenotypes where in_taxon is null or in_taxon = ''").show()

# Check that all phenotypes have at least one causal gene.
