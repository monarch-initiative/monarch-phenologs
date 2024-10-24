import numpy
import pickle
from pathos.multiprocessing import ProcessPool as Pool
from scipy.stats import hypergeom, pearsonr

'''
The phenotype-ortholog matrix is calculated in the prior step.
In this step we generate two additional matrices: The distance matrix and the weight matrix.


'''

# Will need to review this code for efficiency when running on many nodes.
# How can it be optimized to fully utilize the computing cluster?

# Load the resource files.
species_dict = pickle.load(open('../../datasets/utils/species_dict.pkl', 'rb'))
phenotype_list = pickle.load(open('../../datasets/intermediate/gene_candidate/phenotype_list.txt', 'rb'))
test_phenotype_list = pickle.load(open('../../datasets/intermediate/gene_candidate/test_phenotype_list.txt', 'rb'))
ortholog_list = pickle.load(open('../../datasets/intermediate/gene_candidate/ortholog_list.txt', 'rb'))
ortholog_phenotype_matrix = numpy.load('../../datasets/intermediate/gene_candidate/ortholog_phenotype_matrix.npy')

# If testing, set phenotype list to the test list:
phenotype_list = test_phenotype_list

total_phenotypes = len(phenotype_list)
print('INFO: Total number of phenotypes: ' + str(total_phenotypes))
total_orthologs = len(ortholog_list)
print('INFO: Total number of orthologs: ' + str(total_orthologs))
# Have the matrix, need to get the sum of the k (10) nearest neighbors, weight by the pearson correlation coefficient.
# Pearson correlation to determine the k nearest neighbors. So I need to calculate the similarity of phenotypes
# for all pair-wise phenotype combinations. So need a similarity score matrix in addition to the weight matrix.
# Use the hypergeometric CDF to provide scores for the weight matrix.

# Here might be an opportunity to limit the number of phenotypes considered for testing.
# Possible to limit the number of phenotypes but have the full gamut of orthologs for coveage?
# Or would it make more sense just to limit the parsed data in step 01?

distance_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
distance_matrix_comparisons = (len(phenotype_list) * len(phenotype_list))
weight_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
print('Will need to process ' + str(distance_matrix_comparisons) + ' matrix comparisons.')

# Takes ~65 seconds to reach this point.
print('INFO: Assembling phenotype matrix coordinates.')

for phenotype_i in phenotype_list:
    phenotype_i_index = phenotype_list.index(phenotype_i)

    print(
        'INFO: Processing phenotype ' + str(phenotype_i_index + 1) + ' out of ' + str(
            len(phenotype_list)) + '.')
    matrix_coordinates = []
    for phenotype_j in phenotype_list:
        phenotype_j_index = phenotype_list.index(phenotype_j)
        matrix_coordinates.append([phenotype_i_index, phenotype_j_index])

        #without multiprocessing
        # phenotype_index_i = matrix_coordinates[0]
        # phenotype_index_j = matrix_coordinates[1]

        ortholog_counter = 0
        ortholog_match = 0

        (coefficient, p_value) = pearsonr(ortholog_phenotype_matrix[phenotype_i_index],
                                          ortholog_phenotype_matrix[phenotype_j_index])
        for x in range(0, (total_orthologs)):
            if ortholog_phenotype_matrix[phenotype_i_index][x] == 1 and \
                    ortholog_phenotype_matrix[phenotype_j_index][x] == 1:
                ortholog_match += 1
            ortholog_counter += 1

        # N = total number of orthologs shared between species
        # n = nummber of orthologs in species A phenotype
        # m = nummber of orthologs in species B phenotype
        # c = number of common orthologs between phenotypes (ortholog matches)

        m = float(numpy.sum(ortholog_phenotype_matrix[phenotype_i_index]))
        n = float(numpy.sum(ortholog_phenotype_matrix[phenotype_j_index]))
        N = float(total_orthologs)
        c = float(ortholog_match)
        hyp_prob = (hypergeom.cdf(c, N, m, n))

        # Add results to the distance matrix and weight matrix
        distance_matrix[phenotype_i_index][phenotype_j_index] = coefficient
        weight_matrix[phenotype_i_index][phenotype_j_index] = hyp_prob

        # return phenotype_index_i, phenotype_index_j, hyp_prob, coefficient




    print('INFO: Done assembling phenotype matrix coordinates.')
    print('INFO: Starting multiprocessing.')
    # results = m.run(nodes, matrix_coordinates)  # self.multiprocess_matrix_comparisons, matrix_coordinates)
    # results = m.multiprocess_matrix_comparisons(matrix_coordinates)
    print('INFO: Processing results for phenotype ' + str(phenotype_i_index + 1) + ' out of ' + str(
        len(phenotype_list)) + '.')
    print('INFO: Done processing results for phenotype ' + str(phenotype_i_index + 1) + ' out of ' + str(
        len(phenotype_list)) + '.')

# Dump all of the files to disk.
# numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
# numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)
print('INFO: Dumping distance matrix to disk.')
numpy.save('../../datasets/intermediate/gene_candidate/distance_matrix.npy', distance_matrix)
# FYI: The human readable matrix file is 3X the size of the .npy file.
# numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
print('INFO: Dumping weight matrix to disk.')
numpy.save('../../datasets/intermediate/gene_candidate/weight_matrix.npy', weight_matrix)
# numpy.savetxt('inter/phenolog_gene_cand/weight_matrix.txt', weight_matrix)
with open('../../datasets/intermediate/gene_candidate/phenotype_list.txt', 'wb') as handle:
    pickle.dump(phenotype_list, handle)
with open('../../datasets/intermediate/gene_candidate/ortholog_list.txt', 'wb') as handle:
    pickle.dump(ortholog_list, handle)

def multiprocess_matrix_comparisons(self, matrix_coordinates):
        """
        This function processes
        :param matrix_coordinates:
        :return:
        """
        phenotype_index_i = matrix_coordinates[0]
        phenotype_index_j = matrix_coordinates[1]

        ortholog_counter = 0
        ortholog_match = 0

        (coefficient, p_value) = pearsonr(ortholog_phenotype_matrix[phenotype_index_i],
                                          ortholog_phenotype_matrix[phenotype_index_j])
        for x in range(0, (total_orthologs)):
            if ortholog_phenotype_matrix[phenotype_index_i][x] == 1 and \
                    ortholog_phenotype_matrix[phenotype_index_j][x] == 1:
                ortholog_match += 1
            ortholog_counter += 1

        # N = total number of orthologs shared between species
        # n = nummber of orthologs in species A phenotype
        # m = nummber of orthologs in species B phenotype
        # c = number of common orthologs between phenotypes (ortholog matches)

        m = float(numpy.sum(ortholog_phenotype_matrix[phenotype_index_i]))
        n = float(numpy.sum(ortholog_phenotype_matrix[phenotype_index_j]))
        N = float(total_orthologs)
        c = float(ortholog_match)
        hyp_prob = (hypergeom.cdf(c, N, m, n))
        return phenotype_index_i, phenotype_index_j, hyp_prob, coefficient









print('All processing completed.')

