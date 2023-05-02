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

# Set up dataset directories
paths = ['../datasets/sources/monarch_kg', '../datasets/intermediate/human', '../datasets/intermediate/mouse',
         '../datasets/intermediate/rat', '../datasets/intermediate/worm', '../datasets/intermediate/zebrafish',
         '../datasets/intermediate/panther', '../datasets/intermediate/random/fdr/fdr_p_value_lists',
         '../datasets/intermediate/random/human', '../datasets/intermediate/random/mouse',
         '../datasets/intermediate/random/rat', '../datasets/intermediate/random/worm',
         '../datasets/intermediate/random/zebrafish',
         '../datasets/output/phenologs',
         '../datasets/output/gene_candidates']

for path in paths:
    if not os.path.exists(path):
        os.makedirs(path)

# Fetch Monarch KG
# URL = 'https://storage.googleapis.com/monarch-ingest/latest/monarch-kg.tar.gz'
URL = 'https://data.monarchinitiative.org/monarch-kg-dev/latest/monarch-kg.tar.gz'
filename = '../datasets/sources/monarch_kg/' + URL.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(URL)
    f.write(r.content)

# Unpack Monarch KG
monarch_kg = '../datasets/sources/monarch_kg/monarch-kg.tar.gz'
file = tarfile.open(monarch_kg)
file.extractall('../datasets/sources/monarch_kg')
file.close()
# URL = 'https://storage.googleapis.com/monarch-ingest/latest/monarch-kg.tar.gz'
# filename = '../datasets/sources/hpoa/' + URL.split("/")[-1]
# with open(filename, "wb") as f:
#     r = requests.get(URL)
#     f.write(r.content)

print('Setup and data acuisition complete.')