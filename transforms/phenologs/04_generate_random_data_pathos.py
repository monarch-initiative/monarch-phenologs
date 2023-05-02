import random

import pandas as pd
from pathos.multiprocessing import ProcessPool as Pool
import pickle
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)


class myClass:
    def __init__(self):
        pass

    def generate_random_data(self, limit):
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
        -Is the pandas dataframe the best data structure for this process? Would a dataframe with an embedded series be better (one row per phenotype (first column) in dataframe, with genes/orthologs as a series in second column)?
            -> Or perhaps it make more sense to maintain the dict of lists format used previously.
        Question: Is it necessary to create randomized data files that are in both directions?
            For example, I’m currently creating a file for Zebrafish that is zebrafish vs mouse,
            but also creating a random file for mouse that is mouse vs zebrafish.
            Is that duplication necessary?
            Or should I only be creating one randomized file for each pairwise species combination?
            Looks like I previously created files for both directions.
        """
        panther_filepath = "../../datasets/intermediate/panther/panther_orthologs.tsv"

        human_dict = {'species_name': 'human', 'gene_prefix': 'HGNC:',
                      'gene_phenotype_filepath': '../datasets/intermediate/human/human_gene_to_phenotype.tsv',
                      'phenotype_to_ortholog_filepath': '../datasets/intermediate/human/human_phenotype_to_ortholog.pkl',
                      'random_filepath': "../datasets/intermediate/random/human/human_vs_"}
        mouse_dict = {'species_name': 'mouse', 'gene_prefix': 'MGI:',
                      'gene_phenotype_filepath': '../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv',
                      'phenotype_to_ortholog_filepath': '../datasets/intermediate/mouse/mouse_phenotype_to_ortholog.pkl',
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

        for species_a in species_dict:
            for species_b in species_dict:
                if species_a == species_b:
                    pass
                else:
                    phenotype_ortholog_file = species_dict[species_a]['phenotype_to_ortholog_filepath']
                    source_species_name = species_dict[species_a]['species_name']
                    source_gene_prefix = species_dict[species_a]['gene_prefix']
                    target_species_name = species_dict[species_b]['species_name']
                    target_gene_prefix = species_dict[species_b]['gene_prefix']
                    output_filepath = species_dict[species_a]['random_filepath'] + species_dict[species_b]['species_name'] + '_'

                    # Load the phenotype_ortholog pickle file.
                    # This file can be used as the schema for creating the randomized phenotype-ortholog files.
                    data = open(phenotype_ortholog_file, 'rb')
                    phenotype_ortholog_dict = pickle.load(data)

                    # Load common orthologs file for the source and target species.
                    common_orthologs_filepath = "../datasets/intermediate/panther/common_orthologs_" + source_species_name + '_vs_' + target_species_name + '.tsv'
                    common_orthologs = pd.read_csv(common_orthologs_filepath, sep='\t', header=0, low_memory=False)

                    # Have organism structures, have common orthologs, now need to replace each
                    # ortholog associated with a phenotype with a random common ortholog without replacement.
                    # organism1_df = organism1_df.sort_values(by='phenotype', axis=0, ignore_index=True)

                    print('Starting randomized dataset ' + str(limit) + ' for ' + source_species_name + ' vs ' +
                          target_species_name + '.')

                    # Here's an approach using the phenotype-ortholog dicts created in transform 03:
                    '''
                    Start with existing phenotype-ortholog dict
                    Create an empty dict for the randomized data
                    For each phenotype in dict
                        Create a new, shuffled ortholog list 
                        Add the phenotype to the dict
                        For each ortholog associated with the phenotype
                            pop an ortholog off the shuffled list and add to the randomized dict
                    
                    '''
                    randomized_phenotype_ortholog_dict = {}

                    for phenotype in phenotype_ortholog_dict:
                        shuffled_orthologs = common_orthologs.ortholog_id.values.tolist()
                        random.shuffle(shuffled_orthologs)
                        randomized_phenotype_ortholog_dict[phenotype] = []
                        for ortholog in phenotype_ortholog_dict[phenotype]:
                            random_ortholog = shuffled_orthologs.pop()
                            randomized_phenotype_ortholog_dict[phenotype].append(random_ortholog)

                    # print(randomized_phenotype_ortholog_dict)
                    output_file = output_filepath + str(limit) + '.pkl'
                    with open(output_file, 'wb') as handle:
                        pickle.dump(randomized_phenotype_ortholog_dict, handle)
                    print('Completed randomized dataset ' + str(limit) + ' for ' + source_species_name + ' vs ' +
                          target_species_name + ' : ' + output_file)
                    del randomized_phenotype_ortholog_dict, phenotype_ortholog_dict, common_orthologs, shuffled_orthologs
        return

    def run(self, limit, nodes):
        # pool = Pool(nodes=nodes) # Use this one to specify nodes
        pool = Pool() # If nodes not specified, will auto-detect
        pool.map(self.generate_random_data, limit)
        return


if __name__ == '__main__':
    m = myClass()
    nodes = 5
    limit = range(1, 11)
    # limit = range(1, 1001)
    m.run(limit, nodes)
    print('Completed all randomized datasets.')