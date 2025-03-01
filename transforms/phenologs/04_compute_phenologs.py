# General imports
import os
import argparse
import multiprocessing as mp

# Custom imports
from phenologs_utils import (PhenologsSpeciesComparison,
                             divide_workload,
                             initiate_phenologs_species_comparison_configs)


def run_final_phenologs_parallel(configs, num_proc: int = 1):
    """
    This will distribute a list of configs to each cpu core requested (i.e. num_proc argument).
    Each config will run a single species vs species phenologs calculation and write results
    """

    # Deal with - and 0 type edge cases, and instantiate all our objects before running computation
    num_proc = max(1, num_proc)
    run_objs = [PhenologsSpeciesComparison.model_validate(config) for config in configs]

    # Evenly (as possible) divide our data into baskets within a basket (list[list,list,list,...])
    if len(run_objs) > 1:
        div_configs = divide_workload(run_objs, num_proc=num_proc)
    else:
        div_configs = [configs]

    # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
    output = mp.Queue()
    pool = mp.Pool(processes=num_proc)
    results = [pool.apply_async(robj.compute_cross_species_phenologs, args=()) for robj in run_objs]
    output = [ p.get() for p in results ]
    print("- Done!")
    return output


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Computes final phenolog calculations between species. Tables \
                                                      are written for each species comparison based on fdr levels.')

        parser.add_argument("-p", "--project_dir", help="Top most project directory", required=True, type=str)
        parser.add_argument("-c", "--cpu_cores", help="Number of cpu cores to use.", required=False, type=int, default=1)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###
    # Basic run command (Human vs. Mouse for 100 trials across 10 cores)
    ###python phenologs_compute_final_phenologs.py -taxon_ids 9606,10090 -c 10 -p path/to/top_level_project_dir/

    # TO DO: Do we want to have this compute Pearson correlation coeffiecient as well?
    # We are computing the "column vector" distances here from the 2013 methods explanation.
    # This computation can serve both as the function to compute k-nearest neighbors AND the weighting funciton.
    # In other words, we can compute both pearson and hyper geometric distances and mix and match methods
    # for combining via naive bayes (and or additive method) from 2010/2013 papers.  
    taxon_ids, comparison_configs = initiate_phenologs_species_comparison_configs(args)
    run_final_phenologs_parallel(comparison_configs, num_proc=args.cpu_cores)
    
    # Run linearly
    #for config in comparison_configs:
    #    print("- Computing phenologs calculations for {} -- {}".format(config["species_a"], config["species_b"]))
    #    PhenologsSpeciesComparison.parse_obj(config).compute_cross_species_phenologs(outdir)
