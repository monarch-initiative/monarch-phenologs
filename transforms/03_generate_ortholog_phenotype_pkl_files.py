import pandas as pd
import pickle
import numpy as np
from numpy import random
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

def build_ortholog_phenotype_data(panther_file, gene_phenotype_file, output_file):
    # The following steps bring in the Panther data, combine it with the gene-phenotype data,
    # and then selects down to a distinct list of phenotypes and their associated orthologs.
    panther_df = pd.read_csv(panther_file, sep='\t', header=0, low_memory=False)
    panther_df = panther_df[["geneA", "ortholog_id"]]
    panther_df = panther_df.rename(columns={'geneA': 'gene'})
    gene_phenotype_df = pd.read_csv(gene_phenotype_file, sep='\t', header=0, low_memory=False)
    phenotype_ortholog_df = pd.merge(gene_phenotype_df, panther_df, on='gene', how='inner')
    # Align the orthologs with the genes.
    phenotype_ortholog_df = phenotype_ortholog_df[["phenotype", "ortholog_id"]]
    phenotype_ortholog_df = phenotype_ortholog_df.drop_duplicates()

    # TODO: Try converting multiple phenotype-ortholog rows into single phenotype-ortholog series rows.
    # Think this will work better as a dict of lists
    '''
    Iterate through the gene to phenotype file.
    For each phenotype:
        If phenotype not in dict:
            Add phenotype to dict
            Create empty list 
            add ortholog to associated list
        else:
            Check if ortholog in list
                If in list
                    pass
                else 
                    Add to list
            
    '''
    # This will generate a phenotype-to-orthologs dictionary in a 1:n manner,
    # with each phenotype being associated with a distinct list of n orthologs.
    phenotype_ortholog_hash = {}
    for _, row in phenotype_ortholog_df.iterrows():
        if row.phenotype not in phenotype_ortholog_hash:
            phenotype_ortholog_hash[row.phenotype] = [row.ortholog_id]
        elif row.ortholog_id not in phenotype_ortholog_hash[row.phenotype]:
            phenotype_ortholog_hash[row.phenotype].append(row.ortholog_id)
        else:
            print('Phenotype and ortholog already present.')
    # print(phenotype_ortholog_hash)
    with open(output_file, 'wb') as handle:
        pickle.dump(phenotype_ortholog_hash, handle)



    '''
    output_df = phenotype_ortholog_df[['phenotype']]
    output_df = output_df.drop_duplicates()
    output_df['orthologs'] = np.empty((len(output_df), 0)).tolist()
    print(output_df)
    for _, row in phenotype_ortholog_df.iterrows():
        # print(output_df.loc[output_df.phenotype == row.phenotype])
        # print(output_df.loc[output_df.phenotype == row.phenotype])
        output_df.loc[output_df.phenotype == row.phenotype, output_df.orthologs].append(row.ortholog_id)
        # if output_df[] row.phenotype in
        # Check output_df for phenotype
            # If phenotype not fine, print error message to console
        # if phenotype is in output_df
            # check associated series for ortholog
            # if ortholog is present
                #pass
            # if ortholog not present
                # add to series for that row/phenotype


        pass
        # print(row.ortholog_id)
        # print(row.phenotype)

    # print(output_df)
    pd.DataFrame(phenotype_ortholog_df).to_csv(output_file, sep="\t", index=False)
    '''
    return

# Get PANTHER data
# Grab all PANTHER data for now, or just for included species?
# Note, I believe thesis code had to piece together panther links from speciesA -> panther, speciesB -> panther, then lookup ortholog matches.
# Monarch KG has direct speciesA gene -> orthologous to -> speciesB gene.
# For reference, panther edges only exist once per edge, meaning that while there will exist a row for:
# geneA -> orthologous to -> geneB,
# there will not be an equivalent row pointing in the other direction: geneB -> orthologous to -> geneA.
# So in order to have a row for every A to B link, we have to duplicate the dataframe, swap columns, and merge
# OR handle the one-directional nature of the ortholog edges in later steps.
panther_orthologs_filepath = "../datasets/intermediate/panther/panther_orthologs.tsv"

# Build mouse phenotype to ortholog
mouse_gene_to_phenotype_filepath = "../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv"
mouse_phenotype_to_ortholog_filepath = "../datasets/intermediate/mouse/mouse_phenotype_to_ortholog.pkl"
build_ortholog_phenotype_data(panther_orthologs_filepath, mouse_gene_to_phenotype_filepath, mouse_phenotype_to_ortholog_filepath)
print('Mouse  phenotype-to-ortholog complete.')

# Build zebrafish phenotype to ortholog
zebrafish_gene_to_phenotype_filepath = "../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
zebrafish_phenotype_to_ortholog_filepath = "../datasets/intermediate/zebrafish/zebrafish_phenotype_to_ortholog.pkl"
build_ortholog_phenotype_data(panther_orthologs_filepath, zebrafish_gene_to_phenotype_filepath, zebrafish_phenotype_to_ortholog_filepath)
print('Zebrafish phenotype-to-ortholog complete.')

# Build rat phenotype to ortholog
rat_gene_to_phenotype_filepath = "../datasets/intermediate/rat/rat_gene_to_phenotype.tsv"
rat_phenotype_to_ortholog_filepath = "../datasets/intermediate/rat/rat_phenotype_to_ortholog.pkl"
build_ortholog_phenotype_data(panther_orthologs_filepath, rat_gene_to_phenotype_filepath, rat_phenotype_to_ortholog_filepath)
print('Rat phenotype-to-ortholog complete.')

# Build worm phenotype to ortholog
worm_gene_to_phenotype_filepath = "../datasets/intermediate/worm/worm_gene_to_phenotype.tsv"
worm_phenotype_to_ortholog_filepath = "../datasets/intermediate/worm/worm_phenotype_to_ortholog.pkl"
build_ortholog_phenotype_data(panther_orthologs_filepath, worm_gene_to_phenotype_filepath, worm_phenotype_to_ortholog_filepath)
print('Worm phenotype-to-ortholog complete.')

# Build human phenotype to ortholog
human_gene_to_phenotype_filepath = "../datasets/intermediate/human/human_gene_to_phenotype.tsv"
human_phenotype_to_ortholog_filepath = "../datasets/intermediate/human/human_phenotype_to_ortholog.pkl"
build_ortholog_phenotype_data(panther_orthologs_filepath, human_gene_to_phenotype_filepath, human_phenotype_to_ortholog_filepath)
print('Human phenotype-to-ortholog complete.')

print('All ortholog extractions complete.')
