'''
Now that the FDR cutoff has been determined, calculate the final phenologs using the FDR cutoff as the
significance cutoff.

Want an output file containing:
-The cross-species phenotype pairing.
-The calculated p-value score associated with that phenotype pair.
-A flag indicating whether the pair are phenologs/significant based on the FDR cutoff.

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
"""
from scipy.stats import hypergeom, pearsonr
import random
from pathos.multiprocessing import ProcessPool as Pool
import pandas as pd
import pickle
import csv
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)






panther_filepath = "../../datasets/intermediate/panther/panther_orthologs.tsv"
fdr_cutoff_file = '../../datasets/intermediate/random/fdr/fdr_cutoff.pkl'
fdr_cutoff = pickle.load(open(fdr_cutoff_file, 'rb'))
print('Using FDR cutoff value: ' + str(fdr_cutoff) + '.')

# Load species dict.
species_dict = pickle.load(open('../../datasets/utils/species_dict.pkl', 'rb'))


phenologs_df = pd.DataFrame(columns=['Phenotype_A', 'Phenotype_B', 'p_value', 'phenolog_flag'])

total_ortholog_matches = 0
total_ortholog_nonmatches = 0
total_hyp_calcs = 0
# phenolog_p_value_list = []
significant_phenolog_count = 0
total_phenotypes_compared = 0
'''
# Use this setup for species list after testing:
species_list = []
species_list_clone = []
for species in species_dict:
    species_list.append(species)
    species_list_clone.append(species)
species_list.sort()
species_list_clone.sort()
'''

species_list = ['human', 'mouse', 'rat', 'worm', 'zebrafish', 'xenopus']
# species_list = ['worm', 'rat']
species_list.sort()
species_list_clone = ['human', 'mouse', 'rat', 'worm', 'zebrafish', 'xenopus']
# species_list_clone = ['worm', 'rat']
species_list_clone.sort()

with open('../../datasets/output/phenologs/all_phenolog_data.tsv', 'w', newline='') as csvfile:
    all_csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
    with open('../../datasets/output/phenologs/significant_phenolog_data.tsv', 'w', newline='') as csvfile:
        sig_csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
        header_row = ('Species A Phenotype ID', 'Species B Phenotype ID', 'p-value', 'Significance')
        all_csvwriter.writerow(header_row)
        sig_csvwriter.writerow(header_row)

        for species_a in species_list:
            for species_b in species_list_clone:
                if species_a == species_b:
                    pass
                else:
                    print("Starting processing of " + species_a + " vs " + species_b + " phenolog calculations.")
                    counter = 0
                    species_a_phenotype_ortholog_file = species_dict[species_a]['phenotype_to_ortholog_filepath']
                    species_b_phenotype_ortholog_file = species_dict[species_b]['phenotype_to_ortholog_filepath']
                    species_a_phenotype_ortholog_dict = pickle.load(open(species_a_phenotype_ortholog_file, 'rb'))
                    species_b_phenotype_ortholog_dict = pickle.load(open(species_b_phenotype_ortholog_file, 'rb'))
                    species_a_phenotype_count = len(species_a_phenotype_ortholog_dict)
                    print(species_a + " phenotype count: " + str(species_a_phenotype_count))
                    species_b_phenotype_count = len(species_b_phenotype_ortholog_dict)
                    print(species_b + " phenotype count: " + str(species_b_phenotype_count))
                    species_cross_product = species_a_phenotype_count * species_b_phenotype_count
                    print(str(species_cross_product) + ' phenolog calculations to perform.')
                    total_phenotypes_compared += species_cross_product
                    species_a_name = species_dict[species_a]['species_name']
                    species_b_name = species_dict[species_b]['species_name']


                    # Load common orthologs file for the source and target species.
                    common_orthologs_filepath = "../../datasets/intermediate/panther/common_orthologs_" + species_a_name + '_vs_' + species_b_name + '.tsv'
                    common_orthologs = pd.read_csv(common_orthologs_filepath, sep='\t', header=0, low_memory=False)
                    shared_ortholog_count = len(common_orthologs)
                    for species_a_phenotype_id in species_a_phenotype_ortholog_dict:
                        # Phenotype for species A
                        species_a_orthologs = species_a_phenotype_ortholog_dict[species_a_phenotype_id]
                        # print(species_a_orthologs)
                        phenotype_a_ortholog_count = len(species_a_orthologs)
                        for species_b_phenotype_id in species_b_phenotype_ortholog_dict:
                            # Phenotype for species B
                            species_b_orthologs = species_b_phenotype_ortholog_dict[species_b_phenotype_id]
                            ortholog_matches = 0
                            ortholog_non_matches = 0
                            phenotype_b_ortholog_count = len(species_b_orthologs)

                            for species_a_ortholog in species_a_orthologs:
                                # Orthologs for species A
                                for species_b_ortholog in species_b_orthologs:
                                    # Orthologs for species B
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
                                # phenolog_p_value_list.append(prb)
                                total_hyp_calcs += 1
                                if prb <= fdr_cutoff:
                                    significance = 'Significant'
                                    significant_phenolog_count += 1

                                else:
                                    significance = 'Not Significant'

                                # new_row = pd.Series(
                                #     {'Phenotype_A': species_a_phenotype_id, 'Phenotype_B': species_b_phenotype_id,
                                #      'p_value': prb, 'phenolog_flag': significance})
                                # phenologs_df = pd.concat([phenologs_df, new_row.to_frame().T], ignore_index=True)

                            else:
                                prb = 1
                                significance = 'No ortholog matches'
                                new_row = (species_a_phenotype_id, species_b_phenotype_id, prb, significance)
                                # new_row = pd.Series(
                                #     {'Phenotype_A': species_a_phenotype_id, 'Phenotype_B': species_b_phenotype_id,
                                #      'p_value': prb, 'phenolog_flag': significance})
                                # phenologs_df = pd.concat([phenologs_df, new_row.to_frame().T], ignore_index=True)
                            new_row = (species_a_phenotype_id, species_b_phenotype_id, prb, significance)
                            all_csvwriter.writerow(new_row)
                            if significance == 'Significant':
                                sig_csvwriter.writerow(new_row)
                            counter += 1
                            if counter % 100000 == 0:
                                print('Completed phenolog calculation ' + str(counter) + ' of ' + str(
                                    species_cross_product) + '.')
                    print("Completed processing of " + species_a + " vs " + species_b + " phenolog calculations.")
                    print('Total Matches so far: ' + str(total_ortholog_matches))
                    print('Total non-matches so far: ' + str(total_ortholog_nonmatches))
                    print('Total phenotype comparisons so far: ' + str(total_phenotypes_compared))
                    print('Total phenolog calculations so far: ' + str(total_hyp_calcs))
                    print('Total significant phenologs so far: ' + str(significant_phenolog_count))
            species_list_clone.remove(species_a)

            print('All phenologs calculated.')
            print('Total Matches: ' + str(total_ortholog_matches))
            print('Total non-matches: ' + str(total_ortholog_nonmatches))
            print('Total phenotype comparisons: ' + str(total_phenotypes_compared))
            print('Total phenolog calculations: ' + str(total_hyp_calcs))
            print('Total significant phenologs: ' + str(significant_phenolog_count))

# Final output files
# full_output_file = '../../datasets/output/phenologs/all_phenolog_data.tsv'
# pd.DataFrame(phenologs_df).to_csv(full_output_file, sep="\t", index=False)
# significant_phenologs = phenologs_df[(phenologs_df["phenolog_flag"] == 'Significant')]
# significant_phenolog_output_file = '../../datasets/output/phenologs/significant_phenolog_data.tsv'
# pd.DataFrame(significant_phenologs).to_csv(significant_phenolog_output_file, sep="\t", index=False)

print('All processing complete.')

