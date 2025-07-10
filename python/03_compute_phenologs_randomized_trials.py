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
                             RandomSpeciesComparison, 
                             divide_workload,
                             initiate_random_species_comparison_configs,
                             bulk_compute_hyper_geom)


def run_comparisons_parallel(config, n_trials: int = 1, num_proc: int = 1):
    """
    This version will distribute a list of trials to each core requested.
    Hyper geometric calculations are bundled together across different trials within a single core, and then computed
    in bulk (so repeat calculations across different trials are collapsed). These pvalues are then mapped back to each 
    trial's set of parameters and the data written.
    """

    # Deal with - and 0 type edge cases 
    num_proc, n_trials = max(1, num_proc), max(1, n_trials)

    # Properly divide our workload to data
    n_trials = [i for i in range(0, n_trials)]

    # Evenly (as possible) divide our trial numbers into baskets within a basket (list[list,list,list,...])
    if len(n_trials) > 1:
        div_trials = divide_workload(n_trials, num_proc=num_proc)
    else:
        div_trials = [n_trials]
    
    # Instantiate all our objects before running computation
    run_objs = [RandomSpeciesComparison.model_validate(config)
                for i in range(0, len(div_trials))]

    # Setup parallel processing overhead
    # TO DO: Might be nice to limit the number or writers at any given time (more complicated though...)
    #maximum_write = 10
    #semaphore = mp.Semaphore(maximum_write)
    output = mp.Queue()
    pool = mp.Pool(processes=num_proc)

    # Kick off jobs via asynchronous processing
    results = [pool.apply_async(robj.run_randomized_comparison_trials, args=(ddd,)) 
               for robj, ddd in zip(run_objs, div_trials)]
    
    # Retrieve results
    output = [p.get() for p in results]
    print("- Done!")
    return output


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Computes randomized comparison trials between species a vs species b')
        parser.add_argument("-p", "--project_dir", help="Top most project directory", required=True, type=str)
        parser.add_argument("-n", "--num_trials", help="Number of random trials to perform", required=False, type=int, default=1)
        parser.add_argument("-c", "--cpu_cores", help="Number of cpu cores to use.", required=False, type=int, default=1)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    taxon_ids, comparison_configs = initiate_random_species_comparison_configs(args)
    for config in comparison_configs:
        print("- Computing {} random trials between {} -- {}".format(args.num_trials, 
                                                                     config["species_a"], 
                                                                     config["species_b"]))

        run_comparisons_parallel(config, n_trials=args.num_trials, num_proc=args.cpu_cores)