# General imports
import os
import argparse
import requests
import tarfile
import pickle


def download_file_url(url: str, outdir: str, extract_gz: bool = False, overwrite: bool = False):
    """
    Will download file from url to outdir/filename
    filename is generated from the last portion of the url split by "/"
    """
    
    # Download and write file
    filename = os.path.join(outdir, url.split("/")[-1])

    if overwrite != True:
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
                                                      existing data or monarch-kg data if already exists within specified \
                                                      projece_dir')
        parser.add_argument("-p", "--project_dir", help="Directory to write files to", required=True, type=str)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    # KG download nodes, edges paths, phenio
    kg_dir_path = os.path.join(args.project_dir, "monarch_kg")
    kg_edges_path = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")
    phenio_path = os.path.join(args.project_dir, "monarch_kg", "phenio-relation-graph.gz")
    

    # Project top level data directories / structure 
    project_dirs = ["monarch_kg",
                    "species_data",
                    "random_trials", 
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
        URL = 'https://data.monarchinitiative.org/monarch-kg-dev/latest/monarch-kg.tar.gz'
        download_file_url(URL, kg_dir_path, extract_gz=True, overwrite=False)
        print("- Download and upacking of monarch kg succesfull...")
    else:
        print("- Skipping monarch kg download... An edges file already exists at {}".format(kg_edges_path))
    
    # Download phenio relation graph (leave gzip format to save space as we can read through it as is)
    if not os.path.isfile(phenio_path):
        URL = 'https://github.com/monarch-initiative/phenio/releases/latest/download/phenio-relation-graph.gz'
        download_file_url(URL, kg_dir_path, extract_gz=False, overwrite=False)
        print("- Download of phenio relation graph succesfull...")
    else:
        print("- Skipping phenio relation graph download... File already exists at {}".format(phenio_path))