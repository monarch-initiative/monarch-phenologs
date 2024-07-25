import duckdb
from duckdb.typing import *
import os
from scipy.stats import hypergeom, pearsonr

duckdb.sql("CREATE TABLE phenotypes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotypes.csv'") # Load the phenotypes table
duckdb.sql("CREATE TABLE genes AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/genes.csv'") # Load the genes table
duckdb.sql("CREATE TABLE orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/orthologs.csv'") # Load the orthologs table
duckdb.sql("CREATE TABLE phenotype_to_ortholog AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_to_ortholog.csv'") # Load the phenotype-ortholog table

duckdb.sql("CREATE TABLE common_orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_orthologs.csv'") # Load the common orthologs table
duckdb.sql("CREATE TABLE common_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/common_ortholog_counts.csv'") # Load the common ortholog counts table
duckdb.sql("CREATE TABLE total_common_orthologs AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/total_common_orthologs.csv'") # Load the total common ortholog counts table
duckdb.sql("CREATE TABLE phenotype_ortholog_counts AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenotype_ortholog_counts.csv'") # Load the phenotype-ortholog counts table
duckdb.sql("CREATE TABLE phenolog_base AS SELECT * FROM '../../../datasets/intermediate/duckdb_tables/phenolog_base.csv'") # Load the phenologs base table

# For consideration: There may be ways to optimize DuckDB to process data more efficiently e.g. creation of indexes.

# in order to utilize the hypergeometric probabily function, need to create it as a duckDB function:


#prb = float(hypergeom.pmf(c, N, m, n))
def calculate_hypergeometric_probablity(phenotype_a_ortholog_count: int, phenotype_b_ortholog_count: int, shared_ortholog_count: int, ortholog_matches: int):
    m = float(phenotype_b_ortholog_count)
    n = float(phenotype_a_ortholog_count)
    N = float(shared_ortholog_count)
    c = float(ortholog_matches)
    probability = float(hypergeom.pmf(c, N, m, n))
    # print("Phenotype B Ortholog Count: " + str(phenotype_b_ortholog_count) + ", Phenotype A Ortholog Count: " + str(phenotype_a_ortholog_count) + ", Shared ortholog count: "+ str(shared_ortholog_count) + ", Ortholog matches: " + str(ortholog_matches) + ", Probability: " + str(probability))
    return probability

duckdb.create_function("hypergeom_prb", calculate_hypergeometric_probablity, [int, int, int, int], DOUBLE)



# Recommend reviewing whether or not the DOUBLE type is appropriate:
# https://duckdb.org/docs/sql/data_types/numeric
# https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html



print("Phenologs base table: ")
duckdb.sql("SELECT * FROM phenolog_base ").show(max_width=10000, max_rows=10)
# duckdb.sql("UPDATE phenolog_base "
#            "SET weight = 1 ")


duckdb.sql("CREATE TABLE updated_phenolog_base AS "
           "SELECT pb.*, "
           "case when poa.ortholog_count is null then 0 else poa.ortholog_count end as phenotype_a_ortholog_count, "
           "case when pob.ortholog_count is null then 0 else pob.ortholog_count end as phenotype_b_ortholog_count, "
           "case when co.common_ortholog_count is null then 0 else co.common_ortholog_count end as ortholog_matches, "
           "tco.ortholog_count as shared_orthologs, "
           # "null as p_value "
           "FROM phenolog_base pb "
           "LEFT JOIN phenotype_ortholog_counts as poa "
           "ON pb.phenotype_a_id = poa.phenotype_id and pb.phenotype_a_in_taxon = poa.phenotype_taxon_id "
           "LEFT JOIN phenotype_ortholog_counts as pob "
           "ON pb.phenotype_b_id = pob.phenotype_id and pb.phenotype_b_in_taxon = pob.phenotype_taxon_id "
           "LEFT JOIN common_ortholog_counts co "
           "ON pb.phenotype_a_id = co.phenotype_a_id AND pb.phenotype_a_in_taxon = co.phenotype_a_taxon_id "
           "AND pb.phenotype_b_id = co.phenotype_b_id AND pb.phenotype_b_in_taxon = co.phenotype_b_taxon_id "
           "left join total_common_orthologs tco "
           "on pb.phenotype_a_in_taxon = tco.species_a_taxon_id and pb.phenotype_b_in_taxon = tco.species_b_taxon_id "
           )


# phenotype_a_id	phenotype_a_taxon_id	phenotype_b_id	phenotype_b_taxon_id	common_ortholog_count


print("Updated phenologs base table: ")
duckdb.sql("SELECT * FROM updated_phenolog_base ").show(max_width=10000, max_rows=10)

duckdb.sql("ALTER TABLE updated_phenolog_base ADD COLUMN p_value DOUBLE;")

# FOR TESTING: Keep rows that have > 5 ortholog_matches:
# duckdb.sql("DELETE FROM updated_phenolog_base WHERE ortholog_matches < 5;")

print("Updated phenologs base table: ")
duckdb.sql("SELECT * FROM updated_phenolog_base ").show(max_width=10000, max_rows=10)

duckdb.sql("UPDATE updated_phenolog_base "
           "SET p_value = hypergeom_prb(phenotype_b_ortholog_count, phenotype_a_ortholog_count, shared_orthologs, ortholog_matches) ")


print("Updated phenologs base table: populated p-values?")
duckdb.sql("SELECT * FROM updated_phenolog_base ").show(max_width=10000, max_rows=10)

print("Updated phenologs base table: > 0 p-values?")
duckdb.sql("SELECT * FROM updated_phenolog_base where p_value > 0").show(max_width=10000, max_rows=10)


print("Updated phenologs base table, zero ortholog count check: drop from final table?")
duckdb.sql("SELECT * FROM updated_phenolog_base WHERE phenotype_a_ortholog_count = 0 or phenotype_b_ortholog_count = 0").show(max_width=10000, max_rows=10)




print("Updated phenologs base table, common ortholog count check: This should return zero rows.")
duckdb.sql("SELECT * FROM updated_phenolog_base WHERE ortholog_matches > 0 and (phenotype_a_ortholog_count = 0 or phenotype_b_ortholog_count = 0)").show(max_width=10000, max_rows=10)


# save to disk -> Overwrite?
duckdb.sql("COPY updated_phenolog_base TO '../../../datasets/intermediate/duckdb_tables/updated_phenolog_base.csv' (HEADER true, DELIMITER '\t');")


# Xenopus phenotype-ortholog check
print("Xenopus phenotype-ortholog check")
duckdb.sql("SELECT a.phenotype_id as phenotype_a_id, a.phenotype_taxon_id as phenotype_a_taxon_id, a.ortholog_count as phenotype_a_orthologs, b.phenotype_id as phenotype_b_id, b.phenotype_taxon_id as phenotype_b_taxon_id, b.ortholog_count as phenotype_b_orthologs "
           "FROM phenotype_ortholog_counts a "
           "INNER JOIN phenotype_ortholog_counts b "
           "on a.phenotype_id = b.phenotype_id "
           "and a.phenotype_taxon_id = 'NCBITaxon:8355' "
           "and b.phenotype_taxon_id = 'NCBITaxon:8364' "
           "").show(max_width=10000, max_rows=10)


print("Updated phenologs base table: testing")
duckdb.sql("SELECT * FROM updated_phenolog_base "
           "WHERE (phenotype_a_id = 'HP:0005107' or phenotype_b_id = 'HP:0005107') "
           "and (phenotype_a_id = 'XPO:0115436' or phenotype_b_id = 'XPO:0115436')").show(max_width=10000, max_rows=100)

print("For comparison, source phenotype-ortholog counts for each phenotype: ")
duckdb.sql("SELECT * FROM phenotype_ortholog_counts "
           "WHERE phenotype_id in ('HP:0011107', 'XPO:0102085', 'HP:0000920', 'HP:0005107', 'XPO:0115436') ").show(max_width=10000, max_rows=100)

print("For comparison, source common ortholog counts for each phenotype: ")
duckdb.sql("SELECT * FROM common_ortholog_counts "
           "WHERE phenotype_a_id in ('HP:0005107', 'XPO:0115436' ) "
           "and phenotype_b_id in ('HP:0005107', 'XPO:0115436') ").show(max_width=10000, max_rows=100)

print("Common ortholog counts:")
duckdb.sql("SELECT * FROM common_ortholog_counts "
           "WHERE phenotype_a_id in ('HP:0005107', 'XPO:0115436' ) "
           "and phenotype_b_id in ('HP:0005107', 'XPO:0115436') "
           "order by phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id").show(max_width=10000, max_rows=100)

print("Common ortholgs")
duckdb.sql("SELECT * FROM common_orthologs "
           "WHERE phenotype_a_id in ('HP:0005107', 'XPO:0115436' ) "
           "and phenotype_b_id in ('HP:0005107', 'XPO:0115436') "
           "order by phenotype_a_id, phenotype_a_taxon_id, phenotype_b_id, phenotype_b_taxon_id").show(max_width=10000, max_rows=100)

print("Updated phenologs base:")
duckdb.sql("SELECT * FROM updated_phenolog_base "
           "WHERE phenotype_a_id in ('HP:0005107', 'XPO:0115436') "
           "and phenotype_b_id in ('HP:0005107', 'XPO:0115436') "
           "order by phenotype_a_id, phenotype_a_in_taxon, phenotype_b_id, phenotype_b_in_taxon").show(max_width=10000, max_rows=100)

print("For comparison, source phenotype-ortholog counts for each phenotype: ")
duckdb.sql("SELECT * FROM phenotype_ortholog_counts "
           "WHERE phenotype_id in ('HP:0005107', 'XPO:0115436') ").show(max_width=10000, max_rows=100)

print("Checking phenotype ortholgo table for  XPO:0115436")

duckdb.sql("SELECT * FROM phenotype_to_ortholog "
           "WHERE phenotype_id in ('XPO:0115436') ").show(max_width=10000, max_rows=100)