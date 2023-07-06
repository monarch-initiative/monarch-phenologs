import pickle
import numpy

'''
This script generates a sorted list of phenotypes and orthologs to use as indices.
It then generates a phenotype-ortholog matrix, to be used in later calculations.
'''

# Load species dict.
species_dict = pickle.load(open('../../datasets/utils/species_dict.pkl', 'rb'))
print('Assembling ortholog and phenotype lists.')
phenotype_list = []
ortholog_list = []
for species in species_dict:
    species_name = species_dict[species]['species_name']
    print('Starting species ' + species_name + '.')
    # phenotype_ortholog_file = species_dict[species]['phenotype_to_ortholog_filepath']
    phenotype_ortholog_dict = pickle.load(open(species_dict[species]['phenotype_to_ortholog_filepath'], 'rb'))
    for phenotype in phenotype_ortholog_dict:
        if phenotype not in phenotype_list:
            phenotype_list.append(phenotype)
        for ortholog in phenotype_ortholog_dict[phenotype]:
            if ortholog not in ortholog_list:
                ortholog_list.append(ortholog)
    print('Finished species ' + species_name + '.')
    total_phenotypes = len(phenotype_list)
    print('Total number of phenotypes so far: ' + str(total_phenotypes))
    total_orthologs = len(ortholog_list)
    print('Total number of orthologs so far: ' + str(total_orthologs))


phenotype_list.sort()
ortholog_list.sort()
total_phenotypes = len(phenotype_list)
print('INFO: Total number of phenotypes: '+str(total_phenotypes))
total_orthologs = len(ortholog_list)
print('INFO: Total number of orthologs: '+str(total_orthologs))


print('Assembling ortholog-phenotype matrices.')
# Create the ortholog-phenotype matrix as a matrix of zeroes.
ortholog_phenotype_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))
# Now, for each species, enter a 1 at each matrix coordinate corresponding to a phenotype and an associated ortholog.
for species in species_dict:
    species_name = species_dict[species]['species_name']
    phenotype_ortholog_dict = pickle.load(open(species_dict[species]['phenotype_to_ortholog_filepath'], 'rb'))
    print('Adding ' + species_name + ' to the phenotype-ortholog matrix.')
    for phenotype in phenotype_ortholog_dict:
        phenotype_index = phenotype_list.index(phenotype)
        for ortholog in phenotype_ortholog_dict[phenotype]:
            ortholog_index = ortholog_list.index(ortholog)
            ortholog_phenotype_matrix[phenotype_index][ortholog_index] = 1
print('Done assembling phenotype-ortholog matrices.')

# Save all files to disk.
with open('../../datasets/intermediate/gene_candidate/phenotype_list.txt', 'wb') as handle:
    pickle.dump(phenotype_list, handle)
with open('../../datasets/intermediate/gene_candidate/ortholog_list.txt', 'wb') as handle:
    pickle.dump(ortholog_list, handle)
numpy.save('../../datasets/intermediate/gene_candidate/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
numpy.savetxt('../../datasets/intermediate/gene_candidate/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)


print('All processing completed.')