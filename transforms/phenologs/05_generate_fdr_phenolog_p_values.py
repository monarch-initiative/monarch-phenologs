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
panther_filepath = "../../datasets/intermediate/panther/panther_orthologs.tsv"
'''
organism_list = ['human', 'mouse', 'rat', 'worm', 'zebrafish']
'''

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

        # Save p-value list to disc for testing
        p_value_cutoff = phenolog_p_value_list[five_percent_position]
        print(phenolog_p_value_list[1:10])
        print(phenolog_p_value_list[-10:-1])
        print(p_value_cutoff)
        return phenolog_p_value_list

    def stage_phenolog_calculations(self, limit):
        '''
        This function manages running the phenologs calculations for all pair-wise species comparisons using
        randomized datasets for a single

        Possible to just configure by writing an auto-cross-species function?

        :return:
        '''

        # Load species dict.
        species_dict = pickle.load(open('../../datasets/utils/species_dict.pkl', 'rb'))

        # While we needed to create randomized phenotype-ortholog files in each direction for each pair of species,
        # both the processing of these files for the FDR as well as the actual phenologs calculations will only
        # need to be computed once for each species pair. So, the solution here would be to drop out a species from list
        # after each iteration of the 'species_a in species_list loop.
        # Otherwise you are just getting duplicate calculations.
        species_list = []
        species_list_clone = []
        for species in species_dict:
            species_list.append(species)
            species_list_clone.append(species)

        species_list.sort()
        species_list_clone.sort()

        print('Starting species list: ' + str(species_list))
        print('Starting clone species list: ' + str(species_list_clone))
        for species_a in species_list:
            for species_b in species_list_clone:
                if species_a == species_b:
                    pass
                else:
                    phenotype_ortholog_file = species_dict[species_a]['phenotype_to_ortholog_filepath']
                    source_species_name = species_dict[species_a]['species_name']
                    source_gene_prefix = species_dict[species_a]['gene_prefix']
                    # Source and target phenotype-ortholog files should be the random files!
                    source_phenotype_ortholog_file = species_dict[species_a]['phenotype_to_ortholog_filepath']
                    source_random_phenotype_ortholog_file = species_dict[species_a]['random_filepath'] + species_dict[species_b]['species_name'] + '_' + str(limit) + '.pkl'
                    target_species_name = species_dict[species_b]['species_name']
                    target_gene_prefix = species_dict[species_b]['gene_prefix']
                    target_phenotype_ortholog_file = species_dict[species_b]['phenotype_to_ortholog_filepath']
                    target_random_phenotype_ortholog_file = species_dict[species_a]['random_filepath'] + species_dict[species_b]['species_name'] + '_' + str(limit) + '.pkl'
                    output_filepath = species_dict[species_a]['random_filepath'] + species_dict[species_b]['species_name']

                    p_value_list_filepath = "../../datasets/intermediate/random/fdr/fdr_p_value_lists/" + source_species_name + "_vs_" + target_species_name + '_' + str(limit) + ".pkl"
                    # p_value_cutoff_filepath = "../datasets/intermediate/random/fdr/fdr_cutoffs/" + source_species_name + "_vs_" + target_species_name + '_' + limit + ".txt"

                    # Load common orthologs file for the source and target species.
                    common_orthologs_filepath = "../../datasets/intermediate/panther/common_orthologs_" + source_species_name + '_vs_' + target_species_name + '.tsv'
                    common_orthologs = pd.read_csv(common_orthologs_filepath, sep='\t', header=0, low_memory=False)
                    print('Generating phenolog p-value list for ' + source_species_name + ' vs ' + target_species_name + ' ' + str(limit) + '.')
                    phenolog_p_value_list = self.calculate_fdr_from_random_data(source_random_phenotype_ortholog_file, target_random_phenotype_ortholog_file, common_orthologs)
                    # Perhaps writing the phenolog_p_value_list to disk would be appropriate here, should the processing fail?
                    # Do we need the phenolog_p_value_list or just grab the 5% cutoff p-value? Maybe write both?
                    # five_percent_position = round((len(phenolog_p_value_list)) * 0.05)
                    # phenolog_p_value_list.sort(reverse=False)

                    # Save p-value list to disc for testing
                    # p_value_cutoff = phenolog_p_value_list[five_percent_position]
                    with open(p_value_list_filepath, 'wb') as handle:
                        pickle.dump(phenolog_p_value_list, handle)
                    print('Complete generating phenolog p-value list for ' + source_species_name + ' vs ' + target_species_name + ' ' + str(limit) + '.')
                    del phenolog_p_value_list, common_orthologs
            species_list_clone.remove(species_a)
            print('Species list after ' + species_a + 'completed: ' + str(species_list))
            print('Clone species list after ' + species_a + ' completed: ' + str(species_list_clone))
        return


    def run(self, limit, nodes):
        # pool = Pool(nodes=nodes) # Use this one to specify nodes
        pool = Pool() # If nodes not specified, will auto-detect
        pool.map(self.stage_phenolog_calculations, limit)
        return

'''
# Testing old code
zebrafish_file = zebrafish_random_zvm_filepath + '1.pkl'
mouse_file = mouse_random_mvz_filepath + '1.pkl'
source_gene_prefix = mouse_gene_prefix
target_gene_prefix = zebrafish_gene_prefix
# Load orthologs file and select the common orthologs between the source and target species.
# Make sense to go ahead and create all common ortholog files now?
orthologs_df = pd.read_csv(panther_filepath, sep='\t', header=0, low_memory=False)
common_orthologs = orthologs_df[
    (orthologs_df["geneA"].str.contains(source_gene_prefix, regex=True, na=True)) & (
        orthologs_df["geneB"].str.contains(target_gene_prefix, regex=True, na=True))]
common_orthologs = common_orthologs[['ortholog_id']]
common_orthologs = common_orthologs.drop_duplicates()
p_value_list = myClass.calculate_fdr_from_random_data(myClass, mouse_file, zebrafish_file, common_orthologs)
'''

if __name__ == '__main__':
    m = myClass()
    nodes = 5
    limit = range(1, 11)
    # limit = range(1, 2)
    m.run(limit, nodes)
    print('Completed phenologs calculations for randomized datasets.')
