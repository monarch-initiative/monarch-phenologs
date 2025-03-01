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
import argparse
import requests
import tarfile
import pickle

#from phenologs_utils import (species_dict,
#                             phenologs_data_paths)

def download_file_url(url: str, outdir: str, extract_gz: bool = False, overwrite: bool = False):
    """
    Will download file from url to outdir/filename
    filename is generated from the last portion of the url split by "/"
    """
    
    # Download and write file
    filename = os.path.join(outdir, url.split("/")[-1])

    if overwrite == True:
        if os.path.isfile(filename):
            print("- Warning, file {} already exists... Set overwrite to True to download and replace")
            return

    with open(filename, "wb") as f:
        r = requests.get(url)
        f.write(r.content)
    
    # Extract gzip
    if extract_gz != False:
        file = tarfile.open(filename)
        file.extractall(kg_dir_path)
        file.close()


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Performs phenologs setup steps. Makes relevant directory structures, \
                                                      downloads monarch kg, and upacks the kg files. Will not overwrite \
                                                      existing data or monarch-kg data if already exists')
        parser.add_argument("-p", "--project_dir", help="Directory to write files to", required=True, type=str)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    kg_dir_path = os.path.join(args.project_dir, "monarch_kg")
    kg_edges_path = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")
    kg_hp_path = os.path.join(args.project_dir, "monarch_kg", "hp.obo")
    project_dirs = ["monarch_kg",
                    "species_data",
                    "random_trials", 
                    "random_trials_fdr", 
                    "phenologs_results",
                    "leave_one_out_xvalidation"]

    # Create base project directory
    if not os.path.isdir(args.project_dir):
        print("- Creating project directory(ies) at {}".format(args.project_dir))
        os.makedirs(args.project_dir, exist_ok=True)

    # Create dataset directories
    for pdir in project_dirs:
        os.makedirs(os.path.join(args.project_dir, pdir), exist_ok=True)

    # Download and upack monarch-kg
    if not os.path.isfile(kg_edges_path):

        # Fetch Monarch KG and upack (.gz file)
        print("- Downloading and upacking monarch kg to {}".format(kg_dir_path))
        URL = 'https://data.monarchinitiative.org/monarch-kg-dev/latest/monarch-kg.tar.gz'
        download_file_url(URL, kg_dir_path, extract_gz=True, overwrite=False)
        print("- Download and upacking of monarch kg succesfull...")
    else:
        print("- Skipping monarch kg download... An edges file already exists at {}".format(kg_edges_path))
    
    # Human phenotype ontology (This is what the monarch kg uses)
    # Allows for selection of specific phenotype terms based on select parent classes for more granular queries
    if not os.path.isfile(kg_hp_path):

        # Download latest version hp.obo file 
        print("- Downloading HP .obo file to {}".format(kg_hp_path))
        URL = "http://purl.obolibrary.org/obo/hp.obo"
        download_file_url(URL, kg_dir_path, extract_gz=False, overwrite=False)
        print("- Download of hp.obo file succesfull...")
    else:
        print("- Skipping HP ontology download... File already exists at {}".format(kg_hp_path))
    
        


# Note about orthology source(s).. We could use panther orthology connections / tables directly.
# But it seems easier to gather this informaiton from the monarch kg, but will leave link here
# # Fetch Panther data (.gz file)
# URL = 'http://data.pantherdb.org/ftp/generic_mapping/panther_classifications.tar.gz'
# download_file_url(URL, kg_dir_path, extract_gz=True)