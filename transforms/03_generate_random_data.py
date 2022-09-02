import pandas as pd
import numpy
from numpy import random
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)


def generate_random_data(organism1, organism2, df1, df2, orthologs, df_output1, df_output2):
    """
    Inputs:
    -Specify the two species being randomized
    -Take in both species-specific dataframes
    -Take in panther orthologs
    -
    """






    return


"""
This function creates random data sets for the determination of the FDR cutoff for the significant phenologs.
It takes as input the real phenotype-ortholog hash and the common orthologs between two species.
It creates a test phenotype-ortholog hash using the real phenotype-ortholog hash as a guide.
This allows for the creation of a random data set of similar size and shape
(same number of phenotypes with the same number of associated orthologs for each phenotype).

:param pheno_gene_hash: Phenotype-ortholog hash file for a given species.
:param common_orthologs: Table of common orthologs between two species for the random data set.
:param out_dir: Directory for saving the random data set files.
:param limit: Total number of random data sets to create.
:return:


Steps:
Take in both species-specific dataframes
Take in panther orthologs
Create a panther subset that points between the two organisms
Create a randomized dataset for each species where the phenotypes are maintained but genes are replaced with random orthologs


Note: Maintain the shape of the species dataframes, meaning keep the phenotypes but randomly draw the orthologs

*******
How should the orthologs be configured? One row for every ortholog link between pairs of genes, one row for every 
panther family, or multiple rows if more than one gene pair has the same panther family?
Note that thesis work assembled a distinct list of ortholog IDs, so....

"""

mouse_gene_to_phenotype_filepath = "../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv"
mouse_gene_to_phenotype = pd.read_csv(mouse_gene_to_phenotype_filepath, sep='\t', header=0, low_memory=False)
panther_filepath = "../datasets/intermediate/panther/panther_orthologs.tsv"
panther = pd.read_csv(panther_filepath, sep='\t', header=0, low_memory=False)
zebrafish_gene_to_phenotype_filepath = "../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
zebrafish_gene_to_phenotype = pd.read_csv(zebrafish_gene_to_phenotype_filepath, sep='\t', header=0, low_memory=False)

mouse_v_zebrafish_orthologs = panther[(panther["geneA"].str.contains('MGI:', regex=True, na=True)) & (
    panther["geneB"].str.contains('ZFIN:', regex=True, na=True))]

mouse_v_zebrafish_orthologs['ortholog_id'] = mouse_v_zebrafish_orthologs['geneA'] + '_' + mouse_v_zebrafish_orthologs['geneB']
print(mouse_v_zebrafish_orthologs)

# panther_inverse = panther.rename(columns={'subject': 'object', 'object': 'subject'})
# print(panther)
# print(panther_inverse)
# full_panther = pd.concat([panther, panther_inverse])
# full_panther = full_panther.drop_duplicates()