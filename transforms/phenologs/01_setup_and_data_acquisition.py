'''
Purpose: Acquire the source datasets needed for calculating phenologs and gene candidate predictions.

Datasets needed:
-Monarch KG
-Either Mondo->HPO disease to phenotype annotation from Monarch API (or KGX files) or HPOA
-

Rat Genome Database
- mouse Genome Database
- Zebrafish Information Network

Rat
mouse
Zebrafish
Worm
Chicken
Fission yeast?

These data sources can be expanded and the process rerun as additional model organism gene-phenotype annotations + panther orthologs are made available.

'''
import os
import requests
import tarfile
import pickle

# Set up dataset directories
paths = ['../../datasets/sources/monarch_kg', '../../datasets/intermediate/panther',
         '../../datasets/intermediate/human', '../../datasets/intermediate/random/human',
         '../../datasets/intermediate/mouse', '../../datasets/intermediate/random/mouse',
         '../../datasets/intermediate/rat', '../../datasets/intermediate/random/rat',
         '../../datasets/intermediate/worm', '../../datasets/intermediate/random/worm',
         '../../datasets/intermediate/xenopus', '../../datasets/intermediate/random/xenopus',
         '../../datasets/intermediate/zebrafish', '../../datasets/intermediate/random/zebrafish',
         '../../datasets/intermediate/random/fdr/fdr_p_value_lists', '../../datasets/utils',
         '../../datasets/output/phenologs', '../../datasets/output/gene_candidates']

for path in paths:
    if not os.path.exists(path):
        os.makedirs(path)

# Set up and pickle the species dict for use in later scripts.
human_dict = {'species_name': 'human', 'gene_prefix': 'HGNC:', 'phenotype_prefix': 'HP:',
              'gene_phenotype_filepath': '../../datasets/intermediate/human/human_gene_to_phenotype.tsv',
              'phenotype_to_ortholog_filepath': '../../datasets/intermediate/human/human_phenotype_to_ortholog.pkl',
              'random_filepath': "../../datasets/intermediate/random/human/human_vs_"}
mouse_dict = {'species_name': 'mouse', 'gene_prefix': 'MGI:', 'phenotype_prefix': 'MP:',
              'gene_phenotype_filepath': '../../datasets/intermediate/mouse/mouse_gene_to_phenotype.tsv',
              'phenotype_to_ortholog_filepath': '../../datasets/intermediate/mouse/mouse_phenotype_to_ortholog.pkl',
              'random_filepath': "../../datasets/intermediate/random/mouse/mouse_vs_"}
rat_dict = {'species_name': 'rat', 'gene_prefix': 'RGD:', 'phenotype_prefix': 'MP:',
            'gene_phenotype_filepath': '../../datasets/intermediate/rat/rat_gene_to_phenotype.tsv',
            'phenotype_to_ortholog_filepath': '../../datasets/intermediate/rat/rat_phenotype_to_ortholog.pkl',
            'random_filepath': "../../datasets/intermediate/random/rat/rat_vs_"}
worm_dict = {'species_name': 'worm', 'gene_prefix': 'WB:', 'phenotype_prefix': 'WBPhenotype:',
             'gene_phenotype_filepath': '../../datasets/intermediate/worm/worm_gene_to_phenotype.tsv',
             'phenotype_to_ortholog_filepath': '../../datasets/intermediate/worm/worm_phenotype_to_ortholog.pkl',
             'random_filepath': "../../datasets/intermediate/random/worm/worm_vs_"}
xenopus_dict = {'species_name': 'xenopus', 'gene_prefix': 'Xenbase:', 'phenotype_prefix': 'XPO',
                'gene_phenotype_filepath': '../../datasets/intermediate/xenopus/xenopus_gene_to_phenotype.tsv',
                'phenotype_to_ortholog_filepath': '../../datasets/intermediate/xenopus/xenopus_phenotype_to_ortholog.pkl',
                'random_filepath': "../../datasets/intermediate/random/xenopus/xenopus_vs_"}
zebrafish_dict = {'species_name': 'zebrafish', 'gene_prefix': 'ZFIN:', 'phenotype_prefix': 'ZP:',
                  'gene_phenotype_filepath': '../../datasets/intermediate/zebrafish/zebrafish_gene_to_phenotype.tsv',
                  'phenotype_to_ortholog_filepath': '../../datasets/intermediate/zebrafish/zebrafish_phenotype_to_ortholog.pkl',
                  'random_filepath': "../../datasets/intermediate/random/zebrafish/zebrafish_vs_"}
species_dict = {'human': human_dict, 'mouse': mouse_dict, 'rat': rat_dict, 'worm': worm_dict,
                'xenopus': xenopus_dict, 'zebrafish': zebrafish_dict}

species_dict_filepath = '../../datasets/utils/species_dict.pkl'
with open(species_dict_filepath, 'wb') as handle:
    pickle.dump(species_dict, handle)



# Fetch Monarch KG
# URL = 'https://storage.googleapis.com/monarch-ingest/latest/monarch-kg.tar.gz'
URL = 'https://data.monarchinitiative.org/monarch-kg-dev/latest/monarch-kg.tar.gz'
filename = '../../datasets/sources/monarch_kg/' + URL.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(URL)
    f.write(r.content)

# Unpack Monarch KG
monarch_kg = '../../datasets/sources/monarch_kg/monarch-kg.tar.gz'
file = tarfile.open(monarch_kg)
file.extractall('../../datasets/sources/monarch_kg')
file.close()
# URL = 'https://storage.googleapis.com/monarch-ingest/latest/monarch-kg.tar.gz'
# filename = '../datasets/sources/hpoa/' + URL.split("/")[-1]
# with open(filename, "wb") as f:
#     r = requests.get(URL)
#     f.write(r.content)

print('Setup and data acquisition complete.')