import pandas as pd
import numpy
from numpy import random
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)


def generate_random_data(gene_phenotype_file, orthologs_file, source_gene_prefix, target_gene_prefix, output_filepath, limit=1000):
    """
    :param gene_phenotype_file: Dataframe of organism 1's gene to phenotype annotations.
    :param organism2_df: Dataframe of organism 2's gene to phenotype annotations.
    :param orthologs_df: Table of common orthologs between two species for the random data set.
    :param output_file: Directory for saving the random dataset files for organism 1.
    :param organism2_output: Directory for saving the random dataset files for organism 2.
    :param limit: Total number of random data sets to create.
    :return:

    Process:
    -Take in both species-specific dataframes
    -Take in panther orthologs
    -


    Thoughts on randomization:
    -Should a gene be chosen at random from the total distinct gene list, or should we maintain gene frequencies?
    -Previously used panther IDs, but can we use the native gene IDs, but just use them from the ortholog subset? Thinking that panther IDs might be necessary.

    """

    organism1_df = pd.read_csv(gene_phenotype_file, sep='\t', header=0, low_memory=False)
    organism1_df = organism1_df[['phenotype']]
    orthologs_df = pd.read_csv(orthologs_file, sep='\t', header=0, low_memory=False)

    common_orthologs = orthologs_df[(orthologs_df["geneA"].str.contains(source_gene_prefix, regex=True, na=True)) & (orthologs_df["geneB"].str.contains(target_gene_prefix, regex=True, na=True))]

    common_orthologs = common_orthologs[['ortholog_id']]
    common_orthologs = common_orthologs.drop_duplicates()

    # Have organism structures, have common orthologs,
    # now need to replace each gene associated with a phenotype with random ortholog without replacement
    organism1_df = organism1_df.sort_values(by='phenotype', axis=0, ignore_index=True)


    starting_limit = limit
    ortholog_index = 0
    print('Gene-Phenotype annotations: ' + str(len(organism1_df)) + '.')
    while limit > 0:
        print('Starting randomized dataset ' + str(limit) + ' of ' + str(starting_limit) + '.')
        phenotype_ortholog_df = pd.DataFrame(columns=['phenotype', 'ortholog'])
        shuffled_orthologs = common_orthologs.sample(frac=1)
        current_phenotype = ''
        for j in range(len(organism1_df)):
            if organism1_df.loc[j, "phenotype"] == current_phenotype:
                ortholog = shuffled_orthologs.iloc[ortholog_index]
                d = pd.DataFrame([[organism1_df.loc[j, "phenotype"], ortholog['ortholog_id']]], columns=['phenotype', 'ortholog'])
                phenotype_ortholog_df = pd.concat([phenotype_ortholog_df, d])
                ortholog_index += 1
            else:
                shuffled_orthologs = common_orthologs.sample(frac=1)
                ortholog_index = 0
                ortholog = shuffled_orthologs.iloc[ortholog_index]
                d = pd.DataFrame([[organism1_df.loc[j, "phenotype"], ortholog['ortholog_id']]], columns=['phenotype', 'ortholog'])
                phenotype_ortholog_df = pd.concat([phenotype_ortholog_df, d])
                ortholog_index += 1
            current_phenotype = organism1_df.loc[j, "phenotype"]
            # print(j)
        output_file = output_filepath + str(limit) + '.tsv'
        pd.DataFrame(phenotype_ortholog_df).to_csv(output_file, sep="\t", index=False)
        print('Completed randomized dataset ' + str(limit) + ' of ' + str(starting_limit) + '.')
        limit -= 1

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

panther_filepath = "../datasets/intermediate/panther/panther_orthologs.tsv"

human_gene_prefix = 'HGNC:'
human_gene_to_phenotype_filepath = "../datasets/intermediate/human/human_gene_to_phenotype.tsv"
human_random_hvm_filepath = "../datasets/intermediate/random/human/human_hvm_"
human_random_hvr_filepath = "../datasets/intermediate/random/human/human_hvr_"
human_random_hvw_filepath = "../datasets/intermediate/random/human/human_hvw_"
human_random_hvz_filepath = "../datasets/intermediate/random/human/human_hvz_"

mouse_gene_prefix = 'MGI:'
mouse_gene_to_phenotype_filepath = "../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv"
mouse_random_mvh_filepath = "../datasets/intermediate/random/mouse/mouse_mvh_"
mouse_random_mvr_filepath = "../datasets/intermediate/random/mouse/mouse_mvr_"
mouse_random_mvw_filepath = "../datasets/intermediate/random/mouse/mouse_mvw_"
mouse_random_mvz_filepath = "../datasets/intermediate/random/mouse/mouse_mvz_"

rat_gene_prefix = 'RGD:'
rat_gene_to_phenotype_filepath = "../datasets/intermediate/rat/rat_gene_to_phenotype.tsv"
rat_random_rvh_filepath = "../datasets/intermediate/random/rat/rat_rvh_"
rat_random_rvm_filepath = "../datasets/intermediate/random/rat/rat_rvm_"
rat_random_rvw_filepath = "../datasets/intermediate/random/rat/rat_rvw_"
rat_random_rvz_filepath = "../datasets/intermediate/random/rat/rat_rvz_"

worm_gene_prefix = 'WB:'
worm_gene_to_phenotype_filepath = "../datasets/intermediate/worm/worm_gene_to_phenotype.tsv"
worm_random_wvh_filepath = "../datasets/intermediate/random/worm/worm_wvh_"
worm_random_wvm_filepath = "../datasets/intermediate/random/worm/worm_wvm_"
worm_random_wvr_filepath = "../datasets/intermediate/random/worm/worm_wvr_"
worm_random_wvz_filepath = "../datasets/intermediate/random/worm/worm_wvz_"

zebrafish_gene_prefix = 'ZFIN:'
zebrafish_gene_to_phenotype_filepath = "../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
zebrafish_random_zvh_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvh_"
zebrafish_random_zvm_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvm_"
zebrafish_random_zvr_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvr_"
zebrafish_random_zvw_filepath = "../datasets/intermediate/random/zebrafish/zebrafish_zvw_"

# Create all the human randomized datasets:

# generate_random_data(human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, mouse_gene_prefix, human_random_hvm_filepath, limit=1000)
# generate_random_data(human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, rat_gene_prefix, human_random_hvr_filepath, limit=1000)
# generate_random_data(human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, worm_gene_prefix, human_random_hvw_filepath, limit=1000)
# generate_random_data(human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, zebrafish_gene_prefix, human_random_hvz_filepath, limit=1000)

# Create all the mouse randomized datasets:
# generate_random_data(mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, human_gene_prefix, mouse_random_mvh_filepath, limit=1000)
# generate_random_data(mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, rat_gene_prefix, mouse_random_mvr_filepath, limit=1000)
# generate_random_data(mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, worm_gene_prefix, mouse_random_mvw_filepath, limit=1000)
# generate_random_data(mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, zebrafish_gene_prefix, mouse_random_mvz_filepath, limit=1000)  # DONE

# Create all the rat randomized datasets:
# generate_random_data(rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, human_gene_prefix, rat_random_rvh_filepath, limit=1000)
# generate_random_data(rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, mouse_gene_prefix, rat_random_rvm_filepath, limit=1000)
# generate_random_data(rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, worm_gene_prefix, rat_random_rvw_filepath, limit=1000)
# generate_random_data(rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, zebrafish_gene_prefix, rat_random_rvz_filepath, limit=1000)

# Create all the worm randomized datasets:
# generate_random_data(worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, human_gene_prefix, worm_random_wvh_filepath, limit=1000)
# generate_random_data(worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, mouse_gene_prefix, worm_random_wvm_filepath, limit=1000)
# generate_random_data(worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, rat_gene_prefix, worm_random_wvr_filepath, limit=1000)
# generate_random_data(worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, zebrafish_gene_prefix, worm_random_wvz_filepath, limit=1000)

# Create all the zebrafish randomized datasets:
# generate_random_data(zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, human_gene_prefix, zebrafish_random_zvh_filepath, limit=1000)
generate_random_data(zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, mouse_gene_prefix, zebrafish_random_zvm_filepath, limit=1000)
generate_random_data(zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, rat_gene_prefix, zebrafish_random_zvr_filepath, limit=1000)
generate_random_data(zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, worm_gene_prefix, zebrafish_random_zvw_filepath, limit=1000)






# panther_inverse = panther.rename(columns={'subject': 'object', 'object': 'subject'})
# print(panther)
# print(panther_inverse)
# full_panther = pd.concat([panther, panther_inverse])
# full_panther = full_panther.drop_duplicates()