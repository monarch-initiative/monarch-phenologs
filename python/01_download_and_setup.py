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
                                                      existing data or monarch-kg data if already exists within specified \
                                                      projece_dir')
        parser.add_argument("-p", "--project_dir", help="Directory to write files to", required=True, type=str)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    # KG download nodes, edges paths
    kg_dir_path = os.path.join(args.project_dir, "monarch_kg")
    kg_edges_path = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")
    
    # Ontology download paths
    kg_hp_path = os.path.join(args.project_dir, "monarch_kg", "hp.obo")
    kg_mondo_path = os.path.join(args.project_dir, "monarch_kg", "mondo.obo")
    kg_ddpheno_path = os.path.join(args.project_dir, "monarch_kg", "ddpheno.obo")

    # Project top level data directories / structure 
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
    
    # Ontology downloads (download the ones we need to filter terms for. Some ontologies only consist of abnormal terms,
    # so we don't have to do anything up or downstream. But certain ontolgies should be filtered for only relevant terms
    onto_abbrvs = ["hp", "mondo", "mp", "ddpheno"] # These are currently the only ones that need adjustments
    for onto_abb in onto_abbrvs:
        kg_onto_path = os.path.join(args.project_dir, "monarch_kg", "{}.obo".format(onto_abb))
        if not os.path.isfile(kg_onto_path):
            # Download latest version hp.obo file 
            print("- Downloading {} .obo file to {}".format(onto_abb, kg_onto_path))
            URL = "http://purl.obolibrary.org/obo/{}.obo".format(onto_abb)
            download_file_url(URL, kg_dir_path, extract_gz=False, overwrite=False)
            print("- Download of {}.obo file succesfull...".format(onto_abb))
        else:
            print("- Skipping {} ontology download... File already exists at {}".format(onto_abb, kg_onto_path))



# Note about orthology source(s).. We could use panther orthology connections / tables directly.
# But it seems easier to gather this informaiton from the monarch kg, but will leave link here
# # Fetch Panther data (.gz file)
# URL = 'http://data.pantherdb.org/ftp/generic_mapping/panther_classifications.tar.gz'
# download_file_url(URL, kg_dir_path, extract_gz=True)