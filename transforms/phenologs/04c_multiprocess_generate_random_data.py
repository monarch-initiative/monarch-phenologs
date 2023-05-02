import multiprocessing as mp
import sys
import pandas as pd

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

def foo(i):
    panther_filepath = "../../datasets/intermediate/panther/panther_orthologs.tsv"
    source_organism = 'mouse'
    target_organism = 'zebrafish'

    match source_organism:
        case 'mouse':
            source_gene_to_phenotype_filepath = "../../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv"
            source_gene_prefix = 'MGI:'
            random_filepath = "../../datasets/intermediate/random/mouse/"
            random_file_prefix = "mouse_mv"
        case 'zebrafish':
            source_gene_to_phenotype_filepath = "../../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
            source_gene_prefix = 'ZFIN:'
            random_filepath = "../../datasets/intermediate/random/zebrafish/"
            random_file_prefix = "zebrafish_zv"
        case _:
            print(str(source_organism) + " is not currently included as a source organism.")
            sys.exit()

    match target_organism:
        case 'mouse':
            target_gene_prefix = 'MGI:'
            random_file_suffix = "m_"
        case 'zebrafish':
            target_gene_prefix = 'ZFIN:'
            random_file_suffix = "z_"
        case _:
            print(str(target_organism) + " is not currently included as a target organism.")
            sys.exit()

    output_filepath = random_filepath + random_file_prefix + random_file_suffix

    gene_to_phenotype_df = pd.read_csv(source_gene_to_phenotype_filepath, sep='\t', header=0, low_memory=False)
    gene_to_phenotype_df = gene_to_phenotype_df[['phenotype']]
    orthologs_df = pd.read_csv(panther_filepath, sep='\t', header=0, low_memory=False)

    common_orthologs = orthologs_df[(orthologs_df["geneA"].str.contains(source_gene_prefix, regex=True, na=True)) & (
        orthologs_df["geneB"].str.contains(target_gene_prefix, regex=True, na=True))]

    common_orthologs = common_orthologs[['ortholog_id']]
    common_orthologs = common_orthologs.drop_duplicates()

    # Have organism structures, have common orthologs,
    # now need to replace each gene associated with a phenotype with random ortholog without replacement
    gene_to_phenotype_df = gene_to_phenotype_df.sort_values(by='phenotype', axis=0, ignore_index=True)

    ortholog_index = 0
    print('Gene-Phenotype annotations: ' + str(len(gene_to_phenotype_df)) + '.')
    print('Starting ' + str(source_organism) + ' vs ' + str(target_organism) + ' randomized dataset ' + str(i) + '.')
    phenotype_ortholog_df = pd.DataFrame(columns=['phenotype', 'ortholog'])
    shuffled_orthologs = common_orthologs.sample(frac=1)
    current_phenotype = ''
    for j in range(len(gene_to_phenotype_df)):
        if gene_to_phenotype_df.loc[j, "phenotype"] == current_phenotype:
            ortholog = shuffled_orthologs.iloc[ortholog_index]
            d = pd.DataFrame([[gene_to_phenotype_df.loc[j, "phenotype"], ortholog['ortholog_id']]],
                             columns=['phenotype', 'ortholog'])
            phenotype_ortholog_df = pd.concat([phenotype_ortholog_df, d])
            ortholog_index += 1
        else:
            shuffled_orthologs = common_orthologs.sample(frac=1)
            ortholog_index = 0
            ortholog = shuffled_orthologs.iloc[ortholog_index]
            d = pd.DataFrame([[gene_to_phenotype_df.loc[j, "phenotype"], ortholog['ortholog_id']]],
                             columns=['phenotype', 'ortholog'])
            phenotype_ortholog_df = pd.concat([phenotype_ortholog_df, d])
            ortholog_index += 1
        current_phenotype = gene_to_phenotype_df.loc[j, "phenotype"]
        print(j)
    output_file = output_filepath + str(i) + '.tsv'
    pd.DataFrame(phenotype_ortholog_df).to_csv(output_file, sep="\t", index=False)
    print('Completed randomized dataset ' + str(i) + '.')


if __name__ == "__main__":
    mp.set_start_method("fork")
    with mp.Pool() as pool:
        pool.map(foo, range(10))

