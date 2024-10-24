import pickle
import numpy
import re
import csv
import heapq

'''
Now tha the distance matrix and weight matrix have been filled out,
the next step is to create/assemble the phenolog gene candidate prediction matrix.


'''

# Load the existing distance matrix and weight matrix:
distance_matrix = numpy.load('../../datasets/intermediate/gene_candidate/distance_matrix.npy')
weight_matrix = numpy.load('../../datasets/intermediate/gene_candidate/weight_matrix.npy')

# Load the ortholog-phenotype matrix
ortholog_phenotype_matrix = numpy.load('../../datasets/intermediate/gene_candidate/ortholog_phenotype_matrix.npy')

species_dict = pickle.load(open('../../datasets/utils/species_dict.pkl', 'rb'))
phenotype_list = pickle.load(open('../../datasets/intermediate/gene_candidate/phenotype_list.txt', 'rb'))
test_phenotype_list = pickle.load(open('../../datasets/intermediate/gene_candidate/test_phenotype_list.txt', 'rb'))
ortholog_list = pickle.load(open('../../datasets/intermediate/gene_candidate/ortholog_list.txt', 'rb'))

'''
-For each phenotype, need to identify the 10 nearest neighbor phenotyeps
-Should already have the distance matrix values for all phenotype-phenotype pairs.
-So for each phenotype need to identify the top-ten phenotypes according to the distance matrix.



'''






# Old Code



# Previously used hash structures to hold phenotype-gene relationships. Replicate? Are they even used here?
with open('inter/hpo/human_pheno_gene_hash.txt', 'rb') as handle:
    human_phenotype_gene_hash = pickle.load(handle)
with open('inter/mgi/mouse_pheno_gene_hash.txt', 'rb') as handle:
    mouse_phenotype_gene_hash = pickle.load(handle)
with open('inter/zfin/zebrafish_pheno_gene_hash.txt', 'rb') as handle:
    zebrafish_phenotype_gene_hash = pickle.load(handle)
phenotype_gene_hash = {}
phenotype_gene_hash.update(human_phenotype_gene_hash)
phenotype_gene_hash.update(mouse_phenotype_gene_hash)
phenotype_gene_hash.update(zebrafish_phenotype_gene_hash)

with open('inter/hpo/human_phenotype_id_to_label_hash.txt', 'rb') as handle:
    human_phenotype_id_to_label_hash = pickle.load(handle)
with open('inter/mgi/mouse_phenotype_id_to_label_hash.txt', 'rb') as handle:
    mouse_phenotype_id_to_label_hash = pickle.load(handle)
with open('inter/zfin/zebrafish_phenotype_id_to_label_hash.txt', 'rb') as handle:
    zebrafish_phenotype_id_to_label_hash = pickle.load(handle)
phenotype_id_to_label_hash = {}
phenotype_id_to_label_hash.update(human_phenotype_id_to_label_hash)
phenotype_id_to_label_hash.update(mouse_phenotype_id_to_label_hash)
phenotype_id_to_label_hash.update(zebrafish_phenotype_id_to_label_hash)
nearest_neighbor_hash = {}

# Were these problematic phenotypes? Yes, these are phenotypes not under phenotypic abnormality.
# TODO: A better approach than this hard coded approach would be to construct the HPO graph from scratch from the Monarch KG,
# identify the HPO terms that are not a descendant of 'phenotypic abnormality, and remove them.
phenotype_exclusion_subset = ['HP:0003581', 'HP:0003831', 'HP:0011463', 'HP:0003577', 'HP:0003829',
                              'HP:0003593', 'HP:0003587', 'HP:0003621', 'HP:0003584', 'HP:0003596',
                              'HP:0003623', 'HP:0003680', 'HP:0003674', 'HP:0003812', 'HP:0003676',
                              'HP:0003677', 'HP:0003828', 'HP:0011462', 'HP:0003743', 'HP:0001472',
                              'HP:0001425', 'HP:0012275', 'HP:0001423', 'HP:0010982', 'HP:0001417',
                              'HP:0001466', 'HP:0003745', 'HP:0003744', 'HP:0001419', 'HP:0010984',
                              'HP:0001426', 'HP:0000007', 'HP:0001450', 'HP:0001475', 'HP:0001470',
                              'HP:0001444', 'HP:0001427', 'HP:0000006', 'HP:0001428', 'HP:0001452',
                              'HP:0001442']

phenotype_ortholog_prediction_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))

# Want to get the 10 nearest neighbors for a given phenotype.
# Pass in phenotype indice
# Get the slice for the phenotype indice.
# Create a clone of the slice, sort, and get the value of the top k entries.
with open('out/phenolog_gene_cand/nearest_neighbor_phenotypes.txt', 'w', newline='\n') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
    for y in range(0, len(phenotype_list)):

        test_phenotype_id = phenotype_list[y]
        print('Test phenotype: ' + test_phenotype_id)
        if re.match('HP:.*', test_phenotype_id):
            phenotype_filter = 'HP'
        elif re.match('MP:.*', test_phenotype_id):
            phenotype_filter = 'MP'
        elif re.match('ZP:.*', test_phenotype_id):
            phenotype_filter = 'ZP'
        else:
            print('ERROR: Unknown phenotype prefix for ' + str(test_phenotype_id) + '.')
            break

        test_distance_slice = distance_matrix[y]
        # TODO: Adjust this code so that all non-organismal phenotypes are dropped from HPO.
        # The following code will set distance values to zero in the test distance slice
        # if the matching phenotype is from the same species as the test phenotype.
        for x in range(0, len(phenotype_list)):
            if x != y:
                match_phenotype_id = phenotype_list[x]
                match_phenotype_prefix = re.sub(':.*', '', match_phenotype_id)
                if phenotype_filter == match_phenotype_prefix:
                    test_distance_slice[x] = -1
                if match_phenotype_id in phenotype_exclusion_subset:
                    test_distance_slice[x] = -1

        # NOTE: If you want to adjust the number of nearest neighbors, change it here.
        intermediate_nearest_neighbors = heapq.nlargest(11, range(len(test_distance_slice)),
                                                        test_distance_slice.take)
        phenotype_id = phenotype_list[y]
        phenotype_label = phenotype_id_to_label_hash[phenotype_id]
        nearest_neighbors = []
        nearest_neighbor_ids = []
        nearest_neighbor_labels = []
        # This will print out the phenotype IDs for the k nearest neighbors
        for z in intermediate_nearest_neighbors:
            if z != y:
                nearest_neighbors.append(z)
                nearest_neighbor_ids.append(phenotype_list[z])
        for i in nearest_neighbor_ids:
            # nearest_neighbor_label = 0
            nearest_neighbor_labels.append(phenotype_id_to_label_hash[i])
        print('Input phenotype: ' + phenotype_list[y])
        print('Nearest neighbor phenotypes: ' + str(nearest_neighbor_ids) + '.')
        print(nearest_neighbors)
        # for i in nearest_neighbor_ids:
        # nearest_neighbor_labels.append(phenotype_id_to_label_hash[i])

        # For nearest neighbor output file: phenotype_id, phenotype_label, nn-phenotype_ids, nn-phenotype_labels

        if phenotype_id not in nearest_neighbor_hash:
            nearest_neighbor_hash[phenotype_id] = nearest_neighbor_ids
        output_row = (phenotype_id, phenotype_label, nearest_neighbor_ids, nearest_neighbor_labels)
        csvwriter.writerow(output_row)

        # Next I need to take those k nearest neighbor phenotypes, and calculate the probability that
        # ortholog i is associated with phenotype j based on those k phenotypes,
        # using the phenotype matrix and weight matrix.
        # Take the sum

        # i in nearest neighbors is the phenotype index
        for i in nearest_neighbors:
            nearby_phenotype = ortholog_phenotype_matrix[i]
            for j in range(0, len(nearby_phenotype)):
                if ortholog_phenotype_matrix[i][j] != 0:
                    phenotype_ortholog_prediction_matrix[y][j] += weight_matrix[i][j] * \
                                                                  ortholog_phenotype_matrix[i][j]

with open('../../datasets/intermediate/gene_candidate/nearest_neighbor_hash.txt', 'wb') as handle:
    pickle.dump(nearest_neighbor_hash, handle)
numpy.save('../../datasets/intermediate/gene_candidate/phenotype_ortholog_prediction_matrix.npy',
           phenotype_ortholog_prediction_matrix)
numpy.savetxt('../../datasets/intermediate/gene_candidate/phenotype_ortholog_prediction_matrix.txt',
              phenotype_ortholog_prediction_matrix)


