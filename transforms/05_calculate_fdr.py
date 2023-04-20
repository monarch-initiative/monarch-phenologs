from scipy.stats import hypergeom, pearsonr
import random

import pandas as pd
from pathos.multiprocessing import ProcessPool as Pool
import pickle
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)


'''
Using the randomized datasets, we need to calculate a false discovery rate that will be used to identify the
"significant" phenologs. Because this involves processing 1000s of randomized datasets, will need to
utilize python's multiprocessing capabilities or seek out other methods for improving performance.

Caveats:
Need to check and see if there's any issues with p-values being so small that they are effectively 0.

From the paper:
"Significant phenologs were identified at a FDR of 0.05 by ranking real and permuted phenologs on the basis of
the associated hypergeometric probabilities and selecting a threshold of probability where the proportion of
permuted phenologs above the cutoff accounted for 5% of the phenologs."


Steps:
-For each species pair, iterate through the randomized datasets.
-Calculate p-values, and assemble into list? dataframe?.
-Once completed with all phenolog calculations, will need to assemble
-Identify which p-value represents the 5% significance cutoff.


NOTE: Phenolog calculation script: Could that be generalizable enough to create a utility script to be called from multiple processing scripts?
NOTE: This portion of the pipeline will likely consume a large amount of time, so potentially necessary to utilize 
      a server that can parallel process enough threads to allow for completion in a reasonable period of time.
'''

# Would it make sense to move all of these directory labels to a separate file to be referenced by individual scripts?
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

organism_list = ['human', 'mouse', 'rat', 'worm', 'zebrafish']


class myClass:
    def __init__(self):
        pass

# Old code for fdr calculation
    def calculate_fdr_from_random_data(self, species_a_phenotype_ortholog_file, species_b_phenotype_ortholog_file, shared_orthologs):
        """
        This function performs the phenolog calculations between
        the phenotypes of two species from the random data sets.
        For each cross-species pair of phenotypes, the associated orthologs are compared for matches.
        If two phenotypes have one or more matching orthologs, the hypergeometric probability is calculated.
        P-values are added to a list, which is returned after all phenotype comparisons are complete.
        :param species_a_phenotype_ortholog_file: The phenotype-ortholog hash for species A.
        :param species_b_phenotype_ortholog_file: The phenotype-ortholog hash for species B.
        :param shared_orthologs: The file containing the orthologs shared between the two compared species.
        :return: List of p-values from the hypergeometric probability calculation.
        """

        # Testing out


        total_ortholog_matches = 0
        total_ortholog_nonmatches = 0

        total_hyp_calcs = 0
        phenolog_p_value_list = []

        # pickle.load(open(species_a_po_hash, 'rb'))
        species_a_phenotype_ortholog_dict = pickle.load(open(species_a_phenotype_ortholog_file, 'rb'))
        species_b_phenotype_ortholog_dict = pickle.load(open(species_b_phenotype_ortholog_file, 'rb'))
        # specia_a_length = len(species_a_phenotype_ortholog_dict)
        # specia_b_length = len(species_b_phenotype_ortholog_dict)
        # print(specia_a_length)
        # print(specia_b_length)
        shared_ortholog_count = len(shared_orthologs)
        # print(shared_ortholog_count)

        # Iterate through the phenotypes for each species,
        # determining the number of ortholog matches between the orthologs associated with each phenotype.
        for i in species_a_phenotype_ortholog_dict:
            # Phenotype for species A
            species_a_phenotype_id = i
            species_a_orthologs = species_a_phenotype_ortholog_dict[i]
            #print(species_a_orthologs)
            phenotype_a_ortholog_count = len(species_a_orthologs)

            for j in species_b_phenotype_ortholog_dict:
                # Phenotype for species B
                species_b_phenotype_id = j
                species_b_orthologs = species_b_phenotype_ortholog_dict[j]
                ortholog_matches = 0
                ortholog_non_matches = 0
                phenotype_b_ortholog_count = len(species_b_orthologs)
                for k in species_a_orthologs:
                    # Orthologs for species A
                    species_a_ortholog = k
                    for l in species_b_orthologs:
                        # Orthologs for species B
                        species_b_ortholog = l
                        if species_a_ortholog == species_b_ortholog:
                            ortholog_matches += 1
                            total_ortholog_matches += 1
                        else:
                            ortholog_non_matches += 1
                            total_ortholog_nonmatches += 1

                if ortholog_matches > 0:
                    # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
                    # N = total number of orthologs shared between species
                    # n = number of orthologs in species A phenotype
                    # m = number of orthologs in species B phenotype
                    # c = number of common orthologs between phenotypes (ortholog matches)
                    m = float(phenotype_b_ortholog_count)
                    n = float(phenotype_a_ortholog_count)
                    N = float(shared_ortholog_count)
                    c = float(ortholog_matches)
                    prb = float(hypergeom.pmf(c, N, m, n))
                    phenolog_p_value_list.append(prb)
                    total_hyp_calcs += 1
        print('Total Matches: '+str(total_ortholog_matches))
        print('Total non-matches: '+str(total_ortholog_nonmatches))
        print('Total phenolog calculations: '+str(total_hyp_calcs))

        # Perhaps writing the phenolog_p_value_list to disk would be appropriate here, should the processing fail?
        # Do we need the phenolog_p_value_list or just grab the 5% cutoff p-value? Maybe write both?
        five_percent_position = round((len(phenolog_p_value_list))*0.05)
        phenolog_p_value_list.sort(reverse=False)
        p_value_cutoff = phenolog_p_value_list[five_percent_position]
        print(phenolog_p_value_list[1:10])
        print(phenolog_p_value_list[-10:-1])
        print(p_value_cutoff)
        return phenolog_p_value_list


# Testing old code
zebrafish_file = zebrafish_random_zvm_filepath + '1.pkl'
mouse_file = mouse_random_mvz_filepath + '1.pkl'
source_gene_prefix = mouse_gene_prefix
target_gene_prefix = zebrafish_gene_prefix
# Load orthologs file and select the common ortholgs between the source and target species.
orthologs_df = pd.read_csv(panther_filepath, sep='\t', header=0, low_memory=False)
common_orthologs = orthologs_df[
    (orthologs_df["geneA"].str.contains(source_gene_prefix, regex=True, na=True)) & (
        orthologs_df["geneB"].str.contains(target_gene_prefix, regex=True, na=True))]
common_orthologs = common_orthologs[['ortholog_id']]
common_orthologs = common_orthologs.drop_duplicates()
p_value_list = myClass.calculate_fdr_from_random_data(myClass, mouse_file, zebrafish_file, common_orthologs)

