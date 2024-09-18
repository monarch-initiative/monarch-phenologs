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


# Load tables as needed
duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Create the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/genes.csv'") # Create the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Create the genes table
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv'") # Create the phenotype-ortholog table

# duckdb.sql("CREATE TABLE common_orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_orthologs.csv'") # Load the common orthologs table
# duckdb.sql("CREATE TABLE common_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_ortholog_counts.csv'") # Load the common ortholog counts table
# duckdb.sql("CREATE TABLE phenotype_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_counts.csv'") # Load the phenotype-ortholog counts table

duckdb.sql("CREATE TABLE knn_phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/knn_phenotypes.csv'") # Create the knn_phenotypes table
duckdb.sql("CREATE TABLE knn_orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/knn_orthologs.csv'") # Create the knn_orthologs table

#TODO: Update to 05_ table version when available.
duckdb.sql("CREATE TABLE phenolog_base AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/05_updated_phenolog_base.csv'") # Load the phenologs base table

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
           "pearson_coefficient, pearson_p_value "
           "FROM phenolog_base ")


# Grab a subset of the table for testing:

duckdb.sql("CREATE TABLE knn_subset AS "
           "SELECT * FROM knn_table "
           "WHERE phenotaxon_a_id in ('MP:0000137_NCBITaxon:10090', 'HP:0010441_NCBITaxon:9606', 'MP:0006277_NCBITaxon:10090', 'HP:0100856_NCBITaxon:9606') "
           "or phenotaxon_b_id in ('MP:0000137_NCBITaxon:10090', 'HP:0010441_NCBITaxon:9606', 'MP:0006277_NCBITaxon:10090', 'HP:0100856_NCBITaxon:9606') ")




# "row_number() OVER (PARTITION BY list_sort(array_value(concat(pa.phenotype_id, '_', pa.in_taxon), concat(pb.phenotype_id, '_', pb.in_taxon)))) as phenolog_row "
"row_number() OVER (PARTITION BY ORDER BY ) as knn_rank  "

print("kNN Table: ")
duckdb.sql("SELECT * FROM knn_subset").show(max_width=10000, max_rows=10)


duckdb.sql("CREATE TABLE knn_rank AS "
           "SELECT *, "
           "rank() OVER (PARTITION BY phenotaxon_a_id ORDER BY pearson_coefficient desc) as knn_rank "
           "FROM knn_subset ")

print("knn Rank Table: ")
duckdb.sql("SELECT * FROM knn_rank").show(max_width=10000, max_rows=10)



print("Top 20 for phenotaxon ID: MP:0000137_NCBITaxon:10090")
duckdb.sql("SELECT * FROM knn_rank where phenotaxon_a_id in ('MP:0000137_NCBITaxon:10090') "
           "and knn_rank <= 20 "
           "order by knn_rank ").show(max_width=10000, max_rows=100)


duckdb.sql("SELECT * FROM knn_rank where phenotaxon_a_id in ('MP:0000137_NCBITaxon:10090', 'HP:0010441_NCBITaxon:9606', 'MP:0006277_NCBITaxon:10090', 'HP:0100856_NCBITaxon:9606')").show(max_width=10000, max_rows=10)

'''
Okay... What now? 
From here would need to filter the table for the top 10/ <=10 rows, and then we need to assemble the gene candidates.
Gene candidates should be the orthologs associated with the cross-species phenotypes
Need to configure the associated scoring as well. Recall that there was an additive method and a... max method...? Was there a formal choice in the paper?




'''
print("kNN rank table: ")
duckdb.sql("SELECT * FROM knn_rank").show(max_width=10000, max_rows=100)
# Note: using the subset may result in small ranks/no ranks greater than 10

# So, do we want to add the data to the phenolog_base table or keep it separate for now?
'''
Okay, so we have the ranking done. Can filter down to <= 10. Now need to assemble an appropriate table, 
and whether that gets added to the phenologs base table or stay separate for now is fine.
Let's keep them separate for now, maye tie them together later.

Have previously created a knn_phenotypes table with a defined set of 10 columns for NN phenotypes, but need it to be more flexible.
Maybe just start with a simpler table.
Need to link the knn_subset table to the ortholog table, and then the genes table.
Ultimate goal is have:
phenotype
"phenolog" phenotype/cross-species phenotype
orthologs
gene candidate in the starting taxon

JOINS: Go from knn_rank to phenotype to orthologs to genes?

Not sure if we need to filter based on common orthologs, or if that will be managed by the joins...

Need to implement a method for scoring. How did we do that previously? Choosing between summing the pearson coefficients? Average? Max?

How about marking gene candidates that aren't currently associated with the phenotype? Novel gene candidates?

'''

duckdb.sql("CREATE TABLE knn_test AS "
           "SELECT * "
           "FROM knn_rank knn "
           "left join phenotype_to_ortholog p2o "
           "on knn.phenotaxon_b_id = p2o.phenotaxon_id "
           "left join orthologs ol "
           "on p2o.ortholog_id = ol.ortholog_id "
           "AND ol.gene_a_taxon_id = knn.phenotype_a_in_taxon "
           "AND ol.gene_b_taxon_id = knn.phenotype_b_in_taxon "
           "WHERE knn.knn_rank <= 10 ")  # NOTE: Filter here or do a prefilter earlier in the process?
# This will likely explode a bit.... But needs to be done


print("kNN test table: ")
duckdb.sql("SELECT * FROM knn_test").show(max_width=10000, max_rows=30)


# From this table, we should be able to group by the starting phenotaxon and... The ortholog id?

#TODO: Don't think this is right. Need to review the paper on how to assemble the prediction scores.

duckdb.sql("CREATE TABLE grouped_knn_orthologs AS "
           "SELECT phenotype_id, phenotype_taxon_id, phenotaxon_id, ortholog_id, gene_a_id, gene_a_taxon_id, "
           "sum(pearson_coefficient) as pearson_sum, "
           "avg(pearson_coefficient) as pearson_avg, "
           "max(pearson_coefficient) as pearson_max "
           "FROM knn_test "
           "GROUP BY phenotype_id, phenotype_taxon_id, phenotaxon_id, ortholog_id, gene_a_id, gene_a_taxon_id "
           )


print("grouped_knn_orthologs table: ")
duckdb.sql("SELECT * FROM grouped_knn_orthologs  where phenotype_id = 'HP:0003041'").show(max_width=10000, max_rows=30)

# This seems more appropriate for an intermediate table, can apply
duckdb.sql("CREATE TABLE grouped_knn_orthologs_no_genes AS "
           "SELECT phenotype_id, phenotype_taxon_id, phenotaxon_id, ortholog_id,  "
           "sum(pearson_coefficient) as pearson_sum, "
           "avg(pearson_coefficient) as pearson_avg, "
           "max(pearson_coefficient) as pearson_max "
           "FROM knn_test "
           "GROUP BY phenotype_id, phenotype_taxon_id, phenotaxon_id, ortholog_id "
           )


print("grouped_knn_orthologs table without genes: ")
duckdb.sql("SELECT * FROM grouped_knn_orthologs_no_genes where phenotype_id = 'HP:0003041'").show(max_width=10000, max_rows=30)

duckdb.sql("CREATE TABLE ranked_knn_orthologs AS "
           "SELECT *, "
           "rank() OVER (PARTITION BY phenotype_id ORDER BY pearson_sum desc) as pearson_sum_rank, "
           "rank() OVER (PARTITION BY phenotype_id ORDER BY pearson_avg desc) as pearson_avg_rank, "
           "rank() OVER (PARTITION BY phenotype_id ORDER BY pearson_max desc) as pearson_max_rank "
           "FROM grouped_knn_orthologs_no_genes ")

print("ranked_knn_orthologs table: ")
duckdb.sql("SELECT * FROM ranked_knn_orthologs where phenotype_id = 'HP:0003041'").show(max_width=10000, max_rows=30)