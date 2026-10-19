# General imports
import os
import re
import json
import gzip
import shutil
import argparse
import requests
import tarfile

# Where Monarch publishes the kg (releases in monarch-kg/, every build, released or not, in monarch-kg-dev/) and phenio
KG_SITE = "https://data.monarchinitiative.org"
PHENIO_RELEASES = "https://github.com/monarch-initiative/phenio/releases"


def download_file_url(url: str, outdir: str, extract_gz: bool = False, overwrite: bool = False):
    """
    Will download file from url to outdir/filename
    filename is generated from the last portion of the url split by "/"
    """
    
    # Download and write file
    filename = os.path.join(outdir, url.split("/")[-1])

    if overwrite != True:
        if os.path.isfile(filename):
            print("- Warning, file {} already exists... Set overwrite to True to download and replace".format(filename))
            return filename

    # Streamed to disk. A bad status is an error, rather than an html error page saved under the filename
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024*1024):
                f.write(chunk)
    
    # Extract tarball
    if extract_gz != False:
        file = tarfile.open(filename)
        file.extractall(outdir)
        file.close()

    return filename


def find_kg(graph_version: str):
    """
    The folder of a monarch kg version and the text of its metadata.yaml. A date Monarch never released
    (e.g. 2026-03-06) is only in monarch-kg-dev/
    """

    for folder in ("monarch-kg", "monarch-kg-dev"):
        url = "{}/{}/{}/".format(KG_SITE, folder, graph_version)
        r = requests.get(url + "metadata.yaml", timeout=60)
        if r.status_code == 200:
            return url, r.text
        if r.status_code != 404:
            r.raise_for_status()
    raise ValueError("- ERROR: No monarch kg {} in {}/monarch-kg/ or {}/monarch-kg-dev/".format(graph_version, KG_SITE, KG_SITE))


def kg_versions(metadata: str):
    """
    (kg version, phenio version), as dates, from a kg's metadata.yaml. Releases since 2026-05-09 write
        version: '2026-10-07'  ...  - id: phenio ... version: v2026-10-07   (phenio nested under kg-phenio)
    and earlier ones
        kg-version: 2026-03-06  ...  data: phenio: v2025-11-14
    """

    kg = re.search(r"^(?:kg-)?version: '?([\d-]+)", metadata, re.M)
    phenio = (re.search(r"^\s+phenio: '?v?([\d-]+)", metadata, re.M) or
              re.search(r"- id: phenio\n(?:[ \t]+.*\n)*?[ \t]+version: '?v?([\d-]+)", metadata))
    if kg is None or phenio is None:
        raise ValueError("- ERROR: Could not read the kg and phenio versions from the kg's metadata.yaml")
    return kg.group(1), phenio.group(1)


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Performs phenologs setup steps. Makes relevant directory structures, \
                                                      downloads monarch kg and the phenio it was built with, and upacks \
                                                      the kg files. The versions are written to monarch_kg/versions.json \
                                                      for the later steps. Will not overwrite existing data or monarch-kg \
                                                      data if already exists within specified projece_dir')
        parser.add_argument("-p", "--project_dir", help="Directory to write files to", required=True, type=str)
        parser.add_argument("-gv", "--graph_version", help="Version of the monarch knowledge graph to download: 'latest' (the newest release) or a date like '2026-03-06'", required=False, type=str, default="latest")
        parser.add_argument("-pv", "--phenio_version", help="Phenio version (a date) to use instead of the one the kg was built with. Must be same date or earlier than the kg", required=False, type=str, default=None)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    # KG download nodes, edges paths, phenio
    kg_dir_path = os.path.join(args.project_dir, "monarch_kg")
    kg_edges_path = os.path.join(args.project_dir, "monarch_kg", "monarch-kg_edges.tsv")
    phenio_path = os.path.join(args.project_dir, "monarch_kg", "phenio-relation-graph.gz")
    versions_path = os.path.join(args.project_dir, "monarch_kg", "versions.json")
    

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


    # Which kg and phenio the project holds. Settled on the first run and written to monarch_kg/versions.json, where the
    # later steps read them. A rerun keeps them, so 'latest' moving on never mixes two kgs in one project
    if os.path.isfile(versions_path):
        with open(versions_path) as fh:
            versions = json.load(fh)
        for given, key in ((args.graph_version, "graph_version"), (args.phenio_version, "phenio_version")):
            if given not in (None, "latest", versions[key]):
                raise ValueError("- ERROR: {} already holds {} {}. Use a new project directory for {}".format(args.project_dir, key, versions[key], given))

    # A kg put in place by hand (the test slices). There is nothing to look up, so its versions must be given
    elif os.path.isfile(kg_edges_path):
        if args.graph_version == "latest" or args.phenio_version is None:
            raise ValueError("- ERROR: {} already holds a kg but no versions.json. Give the -gv and -pv it is from".format(kg_dir_path))
        versions = dict(graph_version=args.graph_version, phenio_version=args.phenio_version, graph_url=None)

    # Looked up: 'latest' is the newest release, and the kg's metadata.yaml names the phenio it was built with
    else:
        graph_version = args.graph_version
        if graph_version == "latest":
            r = requests.get(KG_SITE + "/monarch-kg/latest/metadata.yaml", timeout=60)
            r.raise_for_status()
            graph_version = kg_versions(r.text)[0]
        graph_url, metadata = find_kg(graph_version)
        versions = dict(graph_version=graph_version, phenio_version=args.phenio_version or kg_versions(metadata)[1], graph_url=graph_url)

    if versions["phenio_version"] > versions["graph_version"]:
        raise ValueError("- ERROR: Phenio version {phenio_version} is more recent than monarch kg version {graph_version}. Please specify a phenio version that is the same date or earlier than the monarch kg version.".format(**versions))
    with open(versions_path, "w") as fh:
        json.dump(versions, fh, indent=1)
    print("- Monarch kg {graph_version} and phenio {phenio_version}".format(**versions))

    # Download and upack monarch-kg, from its dated folder (never latest/, which can move on mid download)
    if not os.path.isfile(kg_edges_path):
        URL = (versions["graph_url"] or find_kg(versions["graph_version"])[0]) + "monarch-kg.tar.gz"
        download_file_url(URL, kg_dir_path, extract_gz=True, overwrite=False)
        print("- Download and upacking of monarch kg succesfull...")
    else:
        print("- Skipping monarch kg download... An edges file already exists at {}".format(kg_edges_path))
    
    # Download phenio relation graph (leave gzip format to save space as we can read through it as is).
    # Releases before 2024-12 shipped it as a tarball of the .tsv instead, so that is fetched and gzipped here
    if not os.path.isfile(phenio_path):
        URL = '{}/download/v{}/phenio-relation-graph.gz'.format(PHENIO_RELEASES, versions["phenio_version"])
        try:
            download_file_url(URL, kg_dir_path, extract_gz=False, overwrite=False)
        except requests.HTTPError:
            print("- No phenio-relation-graph.gz on that release, trying the older phenio-relation-graph.tar.gz...")
            tar_path = download_file_url(URL.replace(".gz", ".tar.gz"), kg_dir_path, extract_gz=True, overwrite=False)
            tsv_path = os.path.join(kg_dir_path, "phenio-relation-graph.tsv")
            with open(tsv_path, "rb") as src, gzip.open(phenio_path, "wb", compresslevel=6) as dst:
                shutil.copyfileobj(src, dst)
            os.remove(tsv_path)
            os.remove(tar_path)
        print("- Download of phenio relation graph succesfull...")
    else:
        print("- Skipping phenio relation graph download... File already exists at {}".format(phenio_path))
