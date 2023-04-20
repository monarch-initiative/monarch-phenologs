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

# Set up dataset directories
paths = ['../datasets/sources/monarch_kg', '../datasets/intermediate/human', '../datasets/intermediate/mouse',
         '../datasets/intermediate/rat', '../datasets/intermediate/worm', '../datasets/intermediate/zebrafish',
         '../datasets/intermediate/panther', '../datasets/intermediate/random', '../datasets/output',]

for path in paths:
    if not os.path.exists(path):
        os.makedirs(path)

URL = 'https://storage.googleapis.com/monarch-ingest/latest/monarch-kg.tar.gz'
filename = '../datasets/sources/monarch_kg/' + URL.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(URL)
    f.write(r.content)

URL = 'https://storage.googleapis.com/monarch-ingest/latest/monarch-kg.tar.gz'
filename = '../datasets/sources/hpoa/' + URL.split("/")[-1]
with open(filename, "wb") as f:
    r = requests.get(URL)
    f.write(r.content)
