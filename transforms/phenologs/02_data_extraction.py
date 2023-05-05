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
Unique tags for 'biolink:orthologous_to' predicate:
['HGNC','MGI','NCBIGene','Xenbase','ZFIN','RGD','WB','FB','SGD','dictyBase']

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
Rat Phenotype: MP

Worm gene: WB:WBGene00003001
Worm phenotype: WBPhenotype:0001191

Xenopus gene: Xenbase:XB-GENE-17331825
Xenopus phenotype: XPO:0115384

Yeast gene:
Yeast phenotype:

Chicken gene:
Chicken phenotype:

# Refer to the Phenologs paper regarding removal of redundant phenotype sets, collapsing of diseases, etc, if applicable.


'''


import csv
import os
import pandas as pd


def get_gene_phenotype_edges_from_kg(edges_file, nodes_file, subject, object, predicate, output_file):
    edges = pd.read_csv(edges_file, sep='\t', header=0, low_memory=False)
    edges = edges[(edges["predicate"] == predicate) & (edges["subject"].str.contains(subject, regex=True, na=True)) & (edges["object"].str.contains(object, regex=True, na=True))]
    # edges = edges[edges["predicate"] == predicate & edges["subject"].str.contains(subject, regex=True, na=True) & edges["object"].str.contains(object, regex=True, na=True)]
    edges = edges[["subject","predicate","object"]]
    edges = edges.rename(columns={'subject': 'gene', 'object': 'phenotype'})
    nodes = pd.read_csv(nodes_file, sep='\t', header=0, low_memory=False)
    nodes = nodes[['id', 'name']]
    gene_nodes = nodes.rename(columns={'id': 'gene', 'name': 'gene_name'})
    phenotype_nodes = nodes.rename(columns={'id': 'phenotype', 'name': 'phenotype_name'})
    edges = pd.merge(edges, gene_nodes, on='gene', how='left')
    edges = pd.merge(edges, phenotype_nodes, on='phenotype', how='left')
    edges = edges.drop_duplicates()
    edges = edges[['gene', 'gene_name', 'predicate', 'phenotype', 'phenotype_name']]
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

def get_panther_edges_from_edges_kg(input_file, subject, object, predicate, output_file):
    edges = pd.read_csv(input_file, sep='\t', header=0, low_memory=False)
    edges = edges[(edges["predicate"] == predicate) & (edges["subject"].str.contains(subject, regex=True, na=True)) & (edges["object"].str.contains(object, regex=True, na=True))]
    # edges = edges[edges["predicate"] == predicate & edges["subject"].str.contains(subject, regex=True, na=True) & edges["object"].str.contains(object, regex=True, na=True)]
    edges = edges[["subject","predicate","object", "has_evidence"]]
    # edges = edges["has_evidence"].str.split(':').str[1]
    edges["has_evidence"] = edges["has_evidence"].str.split(':').str[1]
    # Create a duplicated PANTHER file so that we have rows for orthologs pointing in both directions:
    # As in, have one line representing gene A is orthologous to gene B,
    # and another line representing gene B is orthologous to gene A.
    edges_inverse = edges.rename(columns={'subject': 'geneB', 'object': 'geneA', 'has_evidence': 'ortholog_id'})
    edges = edges.rename(columns={'subject': 'geneA', 'object': 'geneB', 'has_evidence': 'ortholog_id'})
    edges = pd.concat([edges, edges_inverse])
    edges = edges.drop_duplicates()
    pd.DataFrame(edges).to_csv(output_file, sep="\t", index=False)
    return

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
kg_edges = '../datasets/sources/monarch_kg/monarch-kg_edges.tsv'
# edges = pd.read_csv(kg_edges, sep='\t', header=0, low_memory=False)

kg_nodes = '../datasets/sources/monarch_kg/monarch-kg_nodes.tsv'
# nodes = pd.read_csv(kg_nodes, sep='\t', header=0, low_memory=False)

has_phenotype = 'biolink:has_phenotype'
has_ortholog = 'biolink:orthologous_to'



# Get PANTHER data
# Grab all PANTHER data for now, or just for included species?
# Note, I believe thesis code had to piece together panther links from speciesA -> panther, speciesB -> panther, then lookup ortholog matches.
# Monarch KG has direct speciesA gene -> orthologous to -> speciesB gene.
# For reference, panther edges only exist once per edge, meaning that while there will exist a row for:
# geneA -> orthologous to -> geneB,
# there will not be an equivalent row pointing in the other direction: geneB -> orthologous to -> geneA.
# So in order to have a row for every A to B link, we have to duplicate the dataframe, swap columns, and merge
# OR handle the one-directional nature of the ortholog edges in later steps.
panther_gene_prefix = '' #
panther_phenotype_prefix = '' #
panther_orthologs_filepath = "../../datasets/intermediate/panther/panther_orthologs.tsv"
get_panther_edges_from_edges_kg(kg_edges, panther_gene_prefix, panther_phenotype_prefix, has_ortholog, panther_orthologs_filepath)
print('PANTHER ortholog extraction complete.')

orthologs_df = pd.read_csv(panther_orthologs_filepath, sep='\t', header=0, low_memory=False)

# Create 'common orthologs' files.

human_dict = {'species_name': 'human', 'gene_prefix': 'HGNC:', 'phenotype_prefix': 'HP:',
              'gene_phenotype_filepath': '../datasets/intermediate/human/human_gene_to_phenotype.tsv',
              'phenotype_to_ortholog_filepath': '../datasets/intermediate/human/human_phenotype_to_ortholog.pkl',
              'random_filepath': "../datasets/intermediate/random/human/human_vs_"}
mouse_dict = {'species_name': 'mouse', 'gene_prefix': 'MGI:', 'phenotype_prefix': 'MP:',
              'gene_phenotype_filepath': '../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv',
              'phenotype_to_ortholog_filepath': '../datasets/intermediate/mouse/mouse_phenotype_to_ortholog.pkl',
              'random_filepath': "../datasets/intermediate/random/mouse/mouse_vs_"}
rat_dict = {'species_name': 'rat', 'gene_prefix': 'RGD:', 'phenotype_prefix': 'MP:',
            'gene_phenotype_filepath': '../datasets/intermediate/rat/rat_gene_to_phenotype.tsv',
            'phenotype_to_ortholog_filepath': '../datasets/intermediate/rat/rat_phenotype_to_ortholog.pkl',
            'random_filepath': "../datasets/intermediate/random/rat/rat_vs_"}
worm_dict = {'species_name': 'worm', 'gene_prefix': 'WB:', 'phenotype_prefix': 'WBPhenotype:',
             'gene_phenotype_filepath': '../datasets/intermediate/worm/worm_gene_to_phenotype.tsv',
             'phenotype_to_ortholog_filepath': '../datasets/intermediate/worm/worm_phenotype_to_ortholog.pkl',
             'random_filepath': "../datasets/intermediate/random/worm/worm_vs_"}
xenopus_dict = {'species_name': 'xenopus', 'gene_prefix': 'Xenbase:', 'phenotype_prefix': 'XPO',
             'gene_phenotype_filepath': '../datasets/intermediate/xenopus/xenopus_gene_to_phenotype.tsv',
             'phenotype_to_ortholog_filepath': '../datasets/intermediate/xenopus/xenopus_phenotype_to_ortholog.pkl',
             'random_filepath': "../datasets/intermediate/random/xenopus/xenopus_vs_"}
zebrafish_dict = {'species_name': 'zebrafish', 'gene_prefix': 'ZFIN:', 'phenotype_prefix': 'ZP:',
                  'gene_phenotype_filepath': '../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv',
                  'phenotype_to_ortholog_filepath': '../datasets/intermediate/zebrafish/zebrafish_phenotype_to_ortholog.pkl',
                  'random_filepath': "../datasets/intermediate/random/zebrafish/zebrafish_vs_"}
species_dict = {'human': human_dict, 'mouse': mouse_dict, 'rat': rat_dict, 'worm': worm_dict,
                'xenopus' : xenopus_dict, 'zebrafish': zebrafish_dict}

for species_a in species_dict:
    for species_b in species_dict:
        if species_a == species_b:
            pass
        else:
            source_species_name = species_dict[species_a]['species_name']
            source_gene_prefix = species_dict[species_a]['gene_prefix']
            target_species_name = species_dict[species_b]['species_name']
            target_gene_prefix = species_dict[species_b]['gene_prefix']
            common_orthologs_filepath = "../datasets/intermediate/panther/common_orthologs_" + source_species_name + '_vs_' + target_species_name + '.tsv'
            common_orthologs = orthologs_df[
                (orthologs_df["geneA"].str.contains(source_gene_prefix, regex=True, na=True)) & (
                    orthologs_df["geneB"].str.contains(target_gene_prefix, regex=True, na=True))]
            common_orthologs = common_orthologs[['ortholog_id']]
            common_orthologs = common_orthologs.drop_duplicates()
            pd.DataFrame(common_orthologs).to_csv(common_orthologs_filepath, sep="\t", index=False)
print('Common orthologs files created.')

# Get human gene to phenotype TODO: Are there phenotypes with the MONDO prefix as well? Doesn't look like it.
human_gene_prefix = 'HGNC:'
human_phenotype_prefix = 'HP:'
human_gene_to_phenotype_filepath = "../../datasets/intermediate/human/human_gene_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, human_gene_prefix, human_phenotype_prefix, has_phenotype, human_gene_to_phenotype_filepath)
print('Human gene-to-phenotype complete.')

# Get human ORPHA disease to phenotype
human_disease_prefix = 'ORPHA:'
human_phenotype_prefix = 'HP:'
human_orpha_disease_to_phenotype_filepath = "../../datasets/intermediate/human/human_orpha_disease_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, human_disease_prefix, human_phenotype_prefix, has_phenotype, human_orpha_disease_to_phenotype_filepath)
print('Human ORPHA-HP disease-to-phenotype complete.')

# Get human MONDO disease to phenotype
human_disease_prefix = 'MONDO:'
human_phenotype_prefix = 'HP:'
human_mondo_disease_to_phenotype_filepath = "../../datasets/intermediate/human/human_mondo_disease_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, human_disease_prefix, human_phenotype_prefix, has_phenotype, human_mondo_disease_to_phenotype_filepath)
print('Human MONDO-HP disease-to-phenotype complete.')

# Get human MONDO disease to MONDO phenotype
human_disease_prefix = 'MONDO:'
human_phenotype_prefix = 'MONDO:'
human_mondo_disease_to_mondo_phenotype_filepath = "../../datasets/intermediate/human/human_mondo_disease_to_mondo_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, human_disease_prefix, human_phenotype_prefix, has_phenotype, human_mondo_disease_to_mondo_phenotype_filepath)
print('Human MONDO-MONDO disease-to-phenotype complete.')

# Merge human disease to phenotype files:
mondo_to_hp = pd.read_csv(human_mondo_disease_to_phenotype_filepath, sep='\t', header=0, low_memory=False)
mondo_to_mondo = pd.read_csv(human_mondo_disease_to_mondo_phenotype_filepath, sep='\t', header=0, low_memory=False)
all_disease_to_phenotype_filepath = "../../datasets/intermediate/human/human_all_disease_to_phenotype.tsv"
all_human_disease_to_phenotype = pd.concat([mondo_to_hp, mondo_to_mondo])
all_human_disease_to_phenotype = all_human_disease_to_phenotype.drop_duplicates()
pd.DataFrame(all_human_disease_to_phenotype).to_csv(all_disease_to_phenotype_filepath, sep="\t", index=False)
print('Human all disease-to-phenotype merge complete.')

# Get mouse gene to phenotype
mouse_gene_prefix = 'MGI:'
mouse_phenotype_prefix = 'MP:'
mouse_gene_to_phenotype_filepath = "../../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, mouse_gene_prefix, mouse_phenotype_prefix, has_phenotype, mouse_gene_to_phenotype_filepath)
print('Mouse gene-to-phenotype complete.')

# Get zebrafish gene to phenotype
zebrafish_gene_prefix = 'ZFIN:'
zebrafish_phenotype_prefix = 'ZP:'
zebrafish_gene_to_phenotype_filepath = "../../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, zebrafish_gene_prefix, zebrafish_phenotype_prefix, has_phenotype, zebrafish_gene_to_phenotype_filepath)
print('Zebrafish gene-to-phenotype complete.')

# Get rat gene to phenotype -> do we currently not have rat phenotypes?
rat_gene_prefix = 'RGD:'
rat_phenotype_prefix = 'MP:'
rat_gene_to_phenotype_filepath = "../../datasets/intermediate/rat/rat_gene_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, rat_gene_prefix, rat_phenotype_prefix, has_phenotype, rat_gene_to_phenotype_filepath)

# Get worm gene to phenotype
worm_gene_prefix = 'WB:'
worm_phenotype_prefix = 'WBPhenotype:'
worm_gene_to_phenotype_filepath = "../../datasets/intermediate/worm/worm_gene_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, worm_gene_prefix, worm_phenotype_prefix, has_phenotype, worm_gene_to_phenotype_filepath)
print('Worm gene-to-phenotype complete.')

# Get xenopus gene to phenotype
xenopus_gene_prefix = 'Xenbase:'
xenopus_phenotype_prefix = 'XPO:'
xenopus_gene_to_phenotype_filepath = "../../datasets/intermediate/xenopus/xenopus_gene_to_phenotype.tsv"
get_gene_phenotype_edges_from_kg(kg_edges, kg_nodes, xenopus_gene_prefix, xenopus_phenotype_prefix, has_phenotype, xenopus_gene_to_phenotype_filepath)
print('Xenopus gene-to-phenotype complete.')


print('All extractions complete.')