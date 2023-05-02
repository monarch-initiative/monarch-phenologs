'''
Now that the FDR cutoff has been determined, calculate the final phenologs using the FDR cutoff as the
significance cutoff.

Want an output file containing:
-The cross-species phentype pairing.
-The calculated p-value score associated with that phenotype pair.
-A flag indicating whether or not the pair are phenologs/significant based on the FDR cutoff.


Note: Need to make a decision regarding the bidirectional indication of phenolog pairs.
Okay to just include phenologs in one direction like:
Phenotype A, Phenotype B, p-value, signficant/phenolog flag
or show bidirectional (although duplicative) like:
Subject | Predicate | Object | p-value
Phenotype A | has_phenolog | Phenotype B | p-value
Phenotype B | has_phenolog | Phenotype A | p-value
'''

"""
This function performs the phenolog calculations between
the phenotypes of two species from the actual data sets.
For each cross-species pair of phenotypes, the associated orthologs are compared for matches.
If two phenotypes have one or more matching orthologs, the hypergeometric probability is calculated.
P-values are added to a list, which is returned after all phenotype comparisons are complete.
:param species_a_phenotype_ortholog_file: The phenotype-ortholog hash for species A.
:param species_b_phenotype_ortholog_file: The phenotype-ortholog hash for species B.
:param shared_orthologs: The file containing the orthologs shared between the two compared species.
:return: List of p-values from the hypergeometric probability calculation.
"""
from scipy.stats import hypergeom, pearsonr
import random

import pandas as pd
import pickle
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)


panther_filepath = "../../datasets/intermediate/panther/panther_orthologs.tsv"
fdr_cutoff_file = "../../datasets/intermediate/random/fdr/"
fdr_cutoff = 0.0072


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


phenologs_df = pd.DataFrame(columns=['Phenotype_A', 'Phenotype_B', 'p_value', 'phenolog_flag'])

total_ortholog_matches = 0
total_ortholog_nonmatches = 0
total_hyp_calcs = 0
phenolog_p_value_list = []
significant_phenolog_count = 0

species_list = ['human', 'mouse', 'rat', 'worm', 'zebrafish']
species_list.sort()
species_list_clone = ['human', 'mouse', 'rat', 'worm', 'zebrafish']
species_list_clone.sort()
for species_a in species_list:
    for species_b in species_list_clone:
        if species_a == species_b:
            pass
        else:
            print("Starting processing of " + species_a + " vs " + species_b + " phenolog calculations.")
            species_a_phenotype_ortholog_file = species_dict[species_a]['phenotype_to_ortholog_filepath']
            species_b_phenotype_ortholog_file = species_dict[species_b]['phenotype_to_ortholog_filepath']
            species_a_phenotype_ortholog_dict = pickle.load(open(species_a_phenotype_ortholog_file, 'rb'))
            species_b_phenotype_ortholog_dict = pickle.load(open(species_b_phenotype_ortholog_file, 'rb'))

            species_a_name = species_dict[species_a]['species_name']
            species_b_name = species_dict[species_b]['species_name']


            # Load common orthologs file for the source and target species.
            common_orthologs_filepath = "../datasets/intermediate/panther/common_orthologs_" + species_a_name + '_vs_' + species_b_name + '.tsv'
            common_orthologs = pd.read_csv(common_orthologs_filepath, sep='\t', header=0, low_memory=False)
            shared_ortholog_count = len(common_orthologs)
            for i in species_a_phenotype_ortholog_dict:
                # Phenotype for species A
                species_a_phenotype_id = i
                species_a_orthologs = species_a_phenotype_ortholog_dict[i]
                # print(species_a_orthologs)
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
                        if prb <= fdr_cutoff:
                            significance = 'Significant'
                            significant_phenolog_count += 1
                        else:
                            significance = 'Not Significant'

                        new_row = pd.Series({'Phenotype_A': species_a_phenotype_id, 'Phenotype_B': species_b_phenotype_id, 'p_value': prb, 'phenolog_flag': significance})
                        pd.concat([phenologs_df, new_row.to_frame().T], ignore_index=True)

                    else:
                        prb = 1
                        significance = 'No ortholog matches'
                        new_row = pd.Series(
                            {'Phenotype_A': species_a_phenotype_id, 'Phenotype_B': species_b_phenotype_id,
                             'p_value': prb, 'phenolog_flag': significance})
                        pd.concat([phenologs_df, new_row.to_frame().T], ignore_index=True)
            print("Completed processing of " + species_a + " vs " + species_b + " phenolog calculations.")
            print('Total Matches so far: ' + str(total_ortholog_matches))
            print('Total non-matches so far: ' + str(total_ortholog_nonmatches))
            print('Total phenolog calculations so far: ' + str(total_hyp_calcs))
            print('Total significant phenologs so far: ' + str(significant_phenolog_count))
    species_list_clone.remove(species_a)

    print('Total Matches: ' + str(total_ortholog_matches))
    print('Total non-matches: ' + str(total_ortholog_nonmatches))
    print('Total phenolog calculations: ' + str(total_hyp_calcs))
    print('Total significant phenologs: ' + str(significant_phenolog_count))

# Final output
full_output_file = '../datasets/output/phenologs/all_phenolog_data.tsv'
pd.DataFrame(phenologs_df).to_csv(full_output_file, sep="\t", index=False)
significant_phenlogs = phenologs_df[(phenologs_df["phenolog_flag"] == 'Significant')]
significant_phenolog_output_file = '../datasets/output/phenologs/significant_phenolog_data.tsv'
pd.DataFrame(significant_phenlogs).to_csv(significant_phenolog_output_file, sep="\t", index=False)


