import numpy
import pickle
from pathos.multiprocessing import ProcessPool as Pool
from scipy.stats import hypergeom, pearsonr

'''


'''

# Will need to review this code for efficiency when running on many nodes.
# How can it be optimized to fully utilize the computing cluster?


class myClass:
    def __init__(self):
        pass

    def run(self, nodes, matrix_coordinate):
        # pool = Pool(nodes=nodes) # Use this one to specify nodes
        pool = Pool() # If nodes not specified, will auto-detect
        results = pool.map(self.multiprocess_matrix_comparisons, matrix_coordinates, total_orthologs)
        return results

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





if __name__ == '__main__':
    m = myClass()
    nodes = 5
    # ortholog_phenotype_matrix = numpy.load('../../datasets/intermediate/gene_candidate/ortholog_phenotype_matrix.npy')  # Is this needed?
    # Set the number of nearest neighbor phenotypes to consider for predictions.
    k = 11

    # Load the resource files.
    species_dict = pickle.load(open('../../datasets/utils/species_dict.pkl', 'rb'))
    phenotype_list = pickle.load(open('../../datasets/intermediate/gene_candidate/phenotype_list.txt', 'rb'))
    ortholog_list = pickle.load(open('../../datasets/intermediate/gene_candidate/ortholog_list.txt', 'rb'))
    ortholog_phenotype_matrix = numpy.load('../../datasets/intermediate/gene_candidate/ortholog_phenotype_matrix.npy')

    total_phenotypes = len(phenotype_list)
    print('INFO: Total number of phenotypes: ' + str(total_phenotypes))
    total_orthologs = len(ortholog_list)
    print('INFO: Total number of orthologs: ' + str(total_orthologs))
    # Have the matrix, need to get the sum of the k (10) nearest neighbors, weight by the pearson correlation coefficient.
    # Pearson correlation to determine the k nearest neighbors. So I need to calculate the similarity of phenotypes
    # for all pair-wise phenotype combinations. So need a similarity score matrix in addition to the weight matrix.
    # Use the hypergeometric CDF to provide scores for the weight matrix.

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
        print('INFO: Done assembling phenotype matrix coordinates.')
        print('INFO: Starting multiprocessing.')
        results = m.run(nodes, matrix_coordinates)  # self.multiprocess_matrix_comparisons, matrix_coordinates)
        # results = m.multiprocess_matrix_comparisons(matrix_coordinates)
        print('INFO: Processing results for phenotype ' + str(phenotype_i_index + 1) + ' out of ' + str(
            len(phenotype_list)) + '.')
        for x in results:
            (phenotype_i_index, phenotype_index_j, hyp_prob, coefficient) = x
            distance_matrix[phenotype_i_index][phenotype_index_j] = coefficient
            weight_matrix[phenotype_i_index][phenotype_index_j] = hyp_prob
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


    print('All processing completed.')

