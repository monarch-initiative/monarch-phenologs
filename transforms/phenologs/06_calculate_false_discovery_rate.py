'''
Have lists of p-values for each pair-wise comparison of phenotypes (in both directions) for every randomized dataset.
Need to assemble the p-value lists for each set of randomized data (so append the p-values for each bidirectional
pairing of species together) and then grab the value that is at the 5% spot.

Then take the average of those 5% values, which would be the FDR.

'''

from scipy.stats import hypergeom, pearsonr
import random

import pandas as pd
from pathos.multiprocessing import ProcessPool as Pool
import pickle
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
from statistics import mean


#
organism_list = ['human', 'mouse', 'rat', 'worm', 'zebrafish']

human_dict = {'species_name': 'human', 'gene_prefix': 'HGNC:',
              'gene_phenotype_filepath': '../datasets/intermediate/human/human_gene_to_phenotype.tsv',
              'phenotype_to_ortholog_filepath': '../datasets/intermediate/human/human_phenotype_to_ortholog.pkl',
              'random_filepath': "../datasets/intermediate/random/human/human_vs_"}
mouse_dict = {'species_name': 'mouse', 'gene_prefix': 'MGI:',
              'gene_phenotype_filepath': '../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv',
              'phenotype_to_ortholog_filepath': '../datasets/intermediate/rat/rat_phenotype_to_ortholog.pkl',
              'random_filepath': "../datasets/intermediate/random/mouse/mouse_vs_"}
rat_dict = {'species_name': 'rat', 'gene_prefix': 'RGD:',
            'gene_phenotype_filepath': '../datasets/intermediate/rat/rat_gene_to_phenotype.tsv',
            'phenotype_to_ortholog_filepath': '../datasets/intermediate/rat/rat_phenotype_to_ortholog.pkl',
            'random_filepath': "../datasets/intermediate/random/rat/rat_vs_"}
worm_dict = {'species_name': 'worm', 'gene_prefix': 'WB:',
             'gene_phenotype_filepath': '../datasets/intermediate/worm/worm_gene_to_phenotype.tsv',
             'phenotype_to_ortholog_filepath': '../datasets/intermediate/worm/worm_phenotype_to_ortholog.pkl',
             'random_filepath': "../datasets/intermediate/random/worm/worm_vs_"}
zebrafish_dict = {'species_name': 'zebrafish', 'gene_prefix': 'ZFIN:',
                  'gene_phenotype_filepath': '../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv',
                  'phenotype_to_ortholog_filepath': '../datasets/intermediate/zebrafish/zebrafish_phenotype_to_ortholog.pkl',
                  'random_filepath': "../datasets/intermediate/random/zebrafish/zebrafish_vs_"}
species_dict = {'human': human_dict, 'mouse': mouse_dict, 'rat': rat_dict, 'worm': worm_dict,
                'zebrafish': zebrafish_dict}
limit = range(1, 11)
# limit = 1
#
fdr_list = []
species_list = ['human', 'mouse', 'rat', 'worm', 'zebrafish']
species_list.sort()
species_list_clone = ['human', 'mouse', 'rat', 'worm', 'zebrafish']
species_list_clone.sort()
for i in limit:

    p_value_list = []
    print('P-value list: ' + str(p_value_list))
    print('Starting species list: ' + str(species_list))
    for species_a in species_list:
        for species_b in species_list_clone:
            if species_a == species_b:
                pass
            else:
                phenotype_ortholog_file = species_dict[species_a]['phenotype_to_ortholog_filepath']
                source_species_name = species_dict[species_a]['species_name']
                # source_gene_prefix = species_dict[species_a]['gene_prefix']
                # source_phenotype_ortholog_file = species_dict[species_a]['phenotype_to_ortholog_filepath']
                target_species_name = species_dict[species_b]['species_name']
                # target_gene_prefix = species_dict[species_b]['gene_prefix']
                # target_phenotype_ortholog_file = species_dict[species_b]['phenotype_to_ortholog_filepath']
                output_filepath = species_dict[species_a]['random_filepath'] + species_dict[species_b]['species_name']

                p_value_list_filepath = "../datasets/intermediate/random/fdr/fdr_p_value_lists/" + source_species_name + "_vs_" + target_species_name + '_' + str(i) + ".pkl"
                # print('Filepath: ' + p_value_list_filepath)
                p_value_list_file = pickle.load(open(p_value_list_filepath, 'rb'))
                p_value_list.extend(p_value_list_file)
                # print(len(p_value_list))
        species_list_clone.remove(species_a)
        print('Species list after ' + species_a + 'completed: ' + str(species_list))
    print('Assembled P-value list: ' + str(p_value_list[1:10]))
    five_percent_position = round((len(p_value_list)) * 0.05)
    p_value_list.sort(reverse=False)
    p_value_cutoff = p_value_list[five_percent_position]
    fdr_list.append(p_value_cutoff)
    print('FDR List: ' + str(fdr_list))
    # print(p_value_list[1:10])
    # print(p_value_list[-10:-1])
    print('FDR Cutoff value for run ' + str(i) + ': ' + str(p_value_cutoff))


fdr_cutoff_value = mean(fdr_list)
print('Final FDR Cutoff value: ' + str(fdr_cutoff_value))


print('All processing complete.')