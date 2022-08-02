'''
Purpose: Extract individual datasets from the Monarch KG

Datasets needed:
Human: disease to phenotype
Human: gene to phenotype
mouse: gene to phenotype
Zebrafish: gene to phenotype
Rat: gene to phenotype
Worm: gene to phenotype
Yeast: gene to phenotype
Chicken: gene to phenotype

PANTHER: gene to ortholog for all above species


Input: Monarch KG
Output: Dataframes for each species with indicated relations


Unique subject tags for has_phenotype relation:
['WB' 'MGI' 'RGD' 'HGNC' 'ORPHA' 'ZFIN' 'MONDO']

Unique object tags for has_phenotype relation:
['WBPhenotype' 'MP' 'HP' 'ZP' 'MONDO']

Tags:
Human disease: MONDO:0021034,
Human gene: HGNC:7877
Human phenotype: HP:0000001

mouse gene: MGI:87888
mouse phenotype: MP:0002572

Zebrafish gene: ZFIN:ZDB-GENE-210324-7
Zebrafish phenotype: ZP:0002478

Rat gene: RGD:1308009
Rat Phenotype: ?

Worm gene: WB:WBGene00003001
Worm phenotype: WBPhenotype:0001191

Yeast gene:
Yeast phenotype:

Chicken gene:
Chicken phenotype:
'''


import csv
import os
import pandas as pd


def get_edges_from_edges_kg(input_file, subject, object, predicate, output_file):
    edges = pd.read_csv(input_file, sep='\t', header=0, low_memory=False)
    edges = edges[(edges["predicate"] == predicate) & (edges["subject"].str.contains(subject, regex=True, na=True)) & (edges["object"].str.contains(object, regex=True, na=True))]
    # edges = edges[edges["predicate"] == predicate & edges["subject"].str.contains(subject, regex=True, na=True) & edges["object"].str.contains(object, regex=True, na=True)]
    edges = edges[["subject","predicate","object"]]
    edges = edges.drop_duplicates()
    pd.DataFrame(edges).to_csv(output_file, sep="\t", index=False)
    return

'''
Method: Given a subject prefix, object prefix, and a biolink identifier, 
extract all edges from the Monarch KG edges file, 
reduce columns to subject/predicate/object, ensure the rows are distinct, 
and then save to a given file name.

Parameters: 
-subject prefix: Used for string matching of the subject prefix.
-object prefix: Used for string matching of the object prefix.
-predicate: 
-input file
-output file
'''


pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
kg_edges = '../datasets/sources/monarch_kg/monarch-kg/monarch-kg_edges.tsv'
edges = pd.read_csv(kg_edges, sep='\t', header=0, low_memory=False)

kg_nodes = '../datasets/sources/monarch_kg/monarch-kg/monarch-kg_nodes.tsv'
nodes = pd.read_csv(kg_nodes, sep='\t', header=0, low_memory=False)

has_phenotype = 'biolink:has_phenotype'
has_ortholog = 'biolink:orthologous_to'

# Get mouse gene to phenotype
mouse_gene_prefix = 'MGI:'
mouse_phenotype_prefix = 'MP:'
mouse_gene_to_phenotype_filepath = "../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv"
get_edges_from_edges_kg(kg_edges, mouse_gene_prefix, mouse_phenotype_prefix, has_phenotype, mouse_gene_to_phenotype_filepath)

# Get zebrafish gene to phenotype
zebrafish_gene_prefix = 'ZFIN:'
zebrafish_phenotype_prefix = 'ZP:'
zebrafish_gene_to_phenotype_filepath = "../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
get_edges_from_edges_kg(kg_edges, zebrafish_gene_prefix, zebrafish_phenotype_prefix, has_phenotype, zebrafish_gene_to_phenotype_filepath)

# Get rat gene to phenotype -> do we currently have rat phenotypes?
# rat_gene_prefix = 'RGD:'
# rat_phenotype_prefix = 'MP:'
# rat_gene_to_phenotype_filepath = "../datasets/intermediate/rat/rat_gene_to_phenotype.tsv"
# get_edges_from_edges_kg(kg_edges, rat_gene_prefix, rat_phenotype_prefix, has_phenotype, rat_gene_to_phenotype_filepath)

# Get worm gene to phenotype
worm_gene_prefix = 'WB:'
worm_phenotype_prefix = 'WBPhenotype:'
worm_gene_to_phenotype_filepath = "../datasets/intermediate/worm/worm_gene_to_phenotype.tsv"
get_edges_from_edges_kg(kg_edges, worm_gene_prefix, worm_phenotype_prefix, has_phenotype, worm_gene_to_phenotype_filepath)

# Get PANTHER data
# Grab all PANTHER data for now, or just for included species?
# Note, I believe thesis code had to piece together panther links from speciesA -> panther, speciesB -> panther, then lookup ortholog matches.
# Monarch KG has direct speciesA gene -> orthologous to -> speciesB gene.
# For reference, panther edges only exist once per edge, meaning that while there will exist a row for:
# geneA -> orthologous to -> geneB,
# there will not be an equivalent row pointing in the other direction: geneB -> orthologous to -> geneA.
panther_gene_prefix = '' #
panther_phenotype_prefix = '' #
panther_orthologs_filepath = "../datasets/intermediate/panther/panther_orthologs.tsv"
get_edges_from_edges_kg(kg_edges, panther_gene_prefix, panther_phenotype_prefix, has_ortholog, panther_orthologs_filepath)