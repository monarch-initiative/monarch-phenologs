import pandas as pd
import numpy
from numpy import random
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

def build_ortholog_phenotype_data(panther_file, gene_phenotype_file, output_file):
    panther_df = pd.read_csv(panther_file, sep='\t', header=0, low_memory=False)
    panther_df = panther_df[["geneA", "ortholog_id"]]
    panther_df = panther_df.rename(columns={'geneA': 'gene'})
    gene_phenotype_df = pd.read_csv(gene_phenotype_file, sep='\t', header=0, low_memory=False)
    phenotype_ortholog_df = pd.merge(gene_phenotype_df, panther_df, on='gene', how='inner')
    # Align the orthologs with the genes.
    phenotype_ortholog_df = phenotype_ortholog_df[["phenotype", "ortholog_id"]]
    phenotype_ortholog_df = phenotype_ortholog_df.drop_duplicates()
    pd.DataFrame(phenotype_ortholog_df).to_csv(output_file, sep="\t", index=False)
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
mouse_phenotype_to_ortholog_filepath = "../datasets/intermediate/mouse/mouse_phenotype_to_ortholog.tsv"
build_ortholog_phenotype_data(panther_orthologs_filepath, mouse_gene_to_phenotype_filepath, mouse_phenotype_to_ortholog_filepath)
print('Mouse  phenotype-to-ortholog complete.')

# Build zebrafish phenotype to ortholog
zebrafish_gene_to_phenotype_filepath = "../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
zebrafish_phenotype_to_ortholog_filepath = "../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype_to_ortholog.tsv"
build_ortholog_phenotype_data(panther_orthologs_filepath, zebrafish_gene_to_phenotype_filepath, zebrafish_phenotype_to_ortholog_filepath)
print('Zebrafish phenotype-to-ortholog complete.')

# Build rat phenotype to ortholog
rat_gene_to_phenotype_filepath = "../datasets/intermediate/rat/rat_gene_to_phenotype.tsv"
rat_phenotype_to_ortholog_filepath = "../datasets/intermediate/rat/rat_gene_to_phenotype_to_ortholog.tsv"
build_ortholog_phenotype_data(panther_orthologs_filepath, rat_gene_to_phenotype_filepath, rat_phenotype_to_ortholog_filepath)
print('Rat phenotype-to-ortholog complete.')

# Build worm phenotype to ortholog
worm_gene_to_phenotype_filepath = "../datasets/intermediate/worm/worm_gene_to_phenotype.tsv"
worm_phenotype_to_ortholog_filepath = "../datasets/intermediate/worm/worm_gene_to_phenotype_to_ortholog.tsv"
build_ortholog_phenotype_data(panther_orthologs_filepath, worm_gene_to_phenotype_filepath, worm_phenotype_to_ortholog_filepath)
print('Worm phenotype-to-ortholog complete.')

# Build human phenotype to ortholog
human_gene_to_phenotype_filepath = "../datasets/intermediate/human/human_gene_to_phenotype.tsv"
human_phenotype_to_ortholog_filepath = "../datasets/intermediate/human/human_gene_to_phenotype_to_ortholog.tsv"
build_ortholog_phenotype_data(panther_orthologs_filepath, human_gene_to_phenotype_filepath, human_phenotype_to_ortholog_filepath)
print('Human phenotype-to-ortholog complete.')

print('All ortholog extractions complete.')
