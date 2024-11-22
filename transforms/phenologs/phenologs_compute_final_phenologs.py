# General imports
import os
import sys
import argparse
import copy
import pandas as pd
import multiprocessing as mp
from IPython.display import display


# Custom imports
from phenologs_utils import (SpeciesComparison,
                             PhenologsSpeciesComparison,
                             divide_workload,
                             initiate_phenologs_species_comparison_configs,
                             bulk_compute_hyper_geom)



if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Computes final phenolog calculations between species. Tables \
                                                      are written for each species comparison based on fdr levels.')

        parser.add_argument("-p", "--project_dir", help="Top most project directory", required=True, type=str)
        parser.add_argument("-c", "--cpu_cores", help="Number of cores to use.", required=False, type=int, default=1)
        parser.add_argument("-taxon_ids", help="Comma separated list of taxon_ids to use", required=False, type=str)
        parser.add_argument("-all", help="Compare all available species to one another", required=False, action='store_true')
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###
    # Basic run command (Human vs. Mouse for 100 trials across 10 cores)
    ###python phenologs_compute_final_phenologs.py -taxon_ids 9606,10090 -c 10 -p path/to/top_level_project_dir/

    outdir = os.path.join(args.project_dir, "phenologs_results")
    fdr_path = os.path.join(args.project_dir,"random_trials_fdr", "fdr_table.tsv")
    taxon_ids, comparison_configs = initiate_phenologs_species_comparison_configs(args)

    for config in comparison_configs:
        print("- Computing phenologs calculations for {} -- {}".format(config["species_a"], config["species_b"]))
        PhenologsSpeciesComparison.parse_obj(config).compute_cross_species_phenologs(outdir)
