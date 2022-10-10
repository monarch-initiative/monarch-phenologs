import pandas as pd
from pathos.multiprocessing import ProcessPool as Pool

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

        """
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
        # parameters tuple for reference: (gene_phenotype_file, orthologs_file, source_gene_prefix, target_gene_prefix, output_filepath)

        all_params = [(human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, mouse_gene_prefix, human_random_hvm_filepath),
                      (human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, rat_gene_prefix, human_random_hvr_filepath),
                      (human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, worm_gene_prefix, human_random_hvw_filepath),
                      (human_gene_to_phenotype_filepath, panther_filepath, human_gene_prefix, zebrafish_gene_prefix, human_random_hvz_filepath),

                      (mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, human_gene_prefix, mouse_random_mvh_filepath),
                      (mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, rat_gene_prefix, mouse_random_mvr_filepath),
                      (mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, worm_gene_prefix, mouse_random_mvw_filepath),
                      (mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, zebrafish_gene_prefix, mouse_random_mvz_filepath),

                      (rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, human_gene_prefix, rat_random_rvh_filepath),
                      (rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, mouse_gene_prefix, rat_random_rvm_filepath),
                      (rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, worm_gene_prefix, rat_random_rvw_filepath),
                      (rat_gene_to_phenotype_filepath, panther_filepath, rat_gene_prefix, zebrafish_gene_prefix, rat_random_rvz_filepath),

                      (worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, human_gene_prefix, worm_random_wvh_filepath),
                      (worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, mouse_gene_prefix, worm_random_wvm_filepath),
                      (worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, rat_gene_prefix, worm_random_wvr_filepath),
                      (worm_gene_to_phenotype_filepath, panther_filepath, worm_gene_prefix, zebrafish_gene_prefix, worm_random_wvz_filepath),

                      (zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, human_gene_prefix, zebrafish_random_zvh_filepath),
                      (zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, mouse_gene_prefix, zebrafish_random_zvm_filepath),
                      (zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, rat_gene_prefix, zebrafish_random_zvr_filepath),
                      (zebrafish_gene_to_phenotype_filepath, panther_filepath, zebrafish_gene_prefix, worm_gene_prefix, zebrafish_random_zvw_filepath)]



        parameters = [(mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, rat_gene_prefix, mouse_random_mvr_filepath),
                      (mouse_gene_to_phenotype_filepath, panther_filepath, mouse_gene_prefix, worm_gene_prefix, mouse_random_mvw_filepath)]
        # gene_phenotype_file=mouse_gene_to_phenotype_filepath, orthologs_file=panther_filepath, source_gene_prefix=mouse_gene_prefix, target_gene_prefix=rat_gene_prefix, output_filepath=mouse_random_mvr_filepat
        # print(parameters[0])
        for x in all_params:
            gene_phenotype_file = x[0]
            orthologs_file = x[1]
            source_gene_prefix = x[2]
            target_gene_prefix = x[3]
            output_filepath = x[4]

            # print(limit)
            organism1_df = pd.read_csv(gene_phenotype_file, sep='\t', header=0, low_memory=False)
            organism1_df = organism1_df[['phenotype']]
            orthologs_df = pd.read_csv(orthologs_file, sep='\t', header=0, low_memory=False)

            common_orthologs = orthologs_df[
                (orthologs_df["geneA"].str.contains(source_gene_prefix, regex=True, na=True)) & (
                    orthologs_df["geneB"].str.contains(target_gene_prefix, regex=True, na=True))]

            common_orthologs = common_orthologs[['ortholog_id']]
            common_orthologs = common_orthologs.drop_duplicates()

            # Have organism structures, have common orthologs,
            # now need to replace each gene associated with a phenotype with random ortholog without replacement
            organism1_df = organism1_df.sort_values(by='phenotype', axis=0, ignore_index=True)

            starting_limit = limit
            ortholog_index = 0
            print('Starting randomized dataset ' + str(limit) + ' for ' + source_gene_prefix + ' vs ' +
                  target_gene_prefix + '.')
            print('Gene-Phenotype annotations for ' + source_gene_prefix + ' ' + str(len(organism1_df)) + '.')


            phenotype_ortholog_df = pd.DataFrame(columns=['phenotype', 'ortholog'])
            shuffled_orthologs = common_orthologs.sample(frac=1)
            current_phenotype = ''
            for j in range(len(organism1_df)):
                if organism1_df.loc[j, "phenotype"] == current_phenotype:
                    ortholog = shuffled_orthologs.iloc[ortholog_index]
                    d = pd.DataFrame([[organism1_df.loc[j, "phenotype"], ortholog['ortholog_id']]],
                                     columns=['phenotype', 'ortholog'])
                    phenotype_ortholog_df = pd.concat([phenotype_ortholog_df, d])
                    ortholog_index += 1
                else:
                    shuffled_orthologs = common_orthologs.sample(frac=1)
                    ortholog_index = 0
                    ortholog = shuffled_orthologs.iloc[ortholog_index]
                    d = pd.DataFrame([[organism1_df.loc[j, "phenotype"], ortholog['ortholog_id']]],
                                     columns=['phenotype', 'ortholog'])
                    phenotype_ortholog_df = pd.concat([phenotype_ortholog_df, d])
                    ortholog_index += 1
                current_phenotype = organism1_df.loc[j, "phenotype"]
                # print(j)
            output_file = output_filepath + str(limit) + '.tsv'
            pd.DataFrame(phenotype_ortholog_df).to_csv(output_file, sep="\t", index=False)
            print('Completed randomized dataset ' + str(limit) + ' for ' + source_gene_prefix + ' vs ' +
                  target_gene_prefix + '.')


        return

    def run(self, limit, nodes):
        pool = Pool(nodes=nodes)
        pool.map(self.generate_random_data, limit)
        return


if __name__ == '__main__':
    m = myClass()
    nodes = 5
    limit = range(1, 11)
    m.run(limit, nodes)
    print('Completed all randomized datasets.')