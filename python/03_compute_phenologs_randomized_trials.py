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


# EXPERIMENTAL 
def run_comparisons_parallel_v2(config, n_trials: int = 1, num_proc: int = 1):
    """
    This version performs 3 parallelization steps
     - fdr hg param calculations
     - hg tests calulcations
     - writing data
    
    The primary advantage to this version is that one can limit the number of cores writing
    data at one time. Because 3 separate pools of workers are created (1 for each step)

    Note, that 1,000 trials for human vs. mouse will take up 4.2Mb per file x 1,000 = 4.2Gb storage
    So this is the theoretical maximum amount of data we could try to write at one time. 
    """

    # Deal with - and 0 type edge cases 
    num_proc, n_trials = max(1, num_proc), max(1, n_trials)

    # Properly divide our workload to data
    # Integers will automatically be cast to list of integers according to trial number
    if type(n_trials) == type(1):
        n_trials = [i for i in range(0, n_trials)]

    # Evenly (as possible) divide our data into baskets within a basket (list[list,list,list,...])
    if len(n_trials) > 1:
        div_trials = divide_workload(n_trials, num_proc=num_proc)
    else:
        div_trials = [n_trials]
    
    # Instantiate all our objects before running computation
    run_objs = [RandomSpeciesComparison.parse_obj(config) for i in range(0, len(div_trials))]

    ############################################################
    ### Randomized trial / hg parameter set parallel compute ###
    output = mp.Queue()
    pool = mp.Pool(processes=num_proc)

    # Kick off jobs via asynchronous processing
    results = [pool.apply_async(robj.run_randomized_comparison_trials_v2, args=(ddd,)) 
               for robj, ddd in zip(run_objs, div_trials)]
    
    # Retrieve results
    output = [ p.get() for p in results ]
    print("- PART1 COMPLETE...")
    

    ########################################
    ### Hyper geometric parallel compute ###
    # Combine previous output set of parameters into a single set to parallel process
    bulk_hg_params = set()
    for v in output:
        bulk_hg_params = bulk_hg_params | v[1]

    div_hg_params = divide_workload(bulk_hg_params, num_proc=num_proc)

    output_part2 = mp.Queue()
    pool = mp.Pool(processes=num_proc)

    # Kick off jobs via asynchronous processing
    results = [pool.apply_async(bulk_compute_hyper_geom, args=(div_hg_params[i],)) for i in range(0, num_proc)]
    
    # Retrieve results
    output_part2 = [ p.get() for p in results ]
    print("- PART2 COMPLETE...")

    ##############################
    ### Parallel write compute ###
    # Combine previous output into single hg param pvalue map
    hg_pvals = {}
    for v in output_part2:
        hg_pvals.update(v)

    # TO DO (Maybe...): There is probably a less memory intensive way to do this... May or may not be worth it
    # Make copies for parrallel compute
    hg_pval_copies = [copy.copy(hg_pvals) for i in range(0, num_proc)]
    div_write_data = divide_workload([v[0] for v in output], num_proc=num_proc)

    # Can limit the number of "writers" here if need be...
    MAX_WRITE = 10
    output_part3 = mp.Queue()
    pool = mp.Pool(processes=MAX_WRITE)

    # Kick off jobs via asynchronous processing
    results = [pool.apply_async(robj.write_trial_data, args=(ddd, hgpv,)) 
               for robj, ddd, hgpv in zip(run_objs, div_write_data, hg_pval_copies)]
    
    # Retrieve results
    output_part3 = [ p.get() for p in results ]
    pool.close()
    print("- PART3 COMPLETE...")

    return output_part3



def initiate_pairwise_comparison_configs(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    """

    # Ensures part1 & part2 of pipeline have been completed
    check_file = os.path.join(input_args.project_dir, "species_data", "species_information_table.tsv")
    check_outdir = os.path.join(input_args.project_dir, "random_trials")
    
    if not os.path.isfile(check_file):
        print("- ERROR, Project species information table doesn't seem to exist. Exiting...")
        sys.exit()
    
    if not os.path.isdir(check_outdir):
        print("- ERROR, Project random_trials directory doesn't seem to exist. Exiting...")
        sys.exit()
    
    # Figure out which species ids / names we have and which ones are relevant
    species_df = pd.read_csv(check_file, sep='\t')
    species_df = species_df[species_df['Total Phenotype Edges'] > 0] # Only want species whith non zero phenotype information
    

    # Pull out taxon_id information (twice as two separate variables)
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon Label"]))}
    org_taxon_ids = copy.copy(ids_to_name)
    
    
    # Format, and add additional keys to make more friendly to input arguments
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon Label"]))}
    ids_to_name.update({sp_id.split(":")[1]:v for sp_id,v in ids_to_name.items()}) # Add keys without NCBITaxon: prefix
    display(species_df)

    # Simple table mapping species name to the total number of genes (with at least one phenotype term associated with them)
    species_g_count = {t:int(v) for t,v in zip(list(species_df["Taxon Label"]), list(species_df["Genes >= 1 Phenotype"]))}

    # Figure out which species ids we are tasked with comparing to one another
    valid_species_ids = set(ids_to_name.keys())
    if not input_args.all:
        if not input_args.taxon_ids:
            print('- ERROR, -all or -taxon_ids argument must be specified. Exiting...')
            sys.exit()
        
        t_ids = input_args.taxon_ids.split(",")
        tot_species = len(t_ids)
        print(valid_species_ids)

        # Ensure input is compatible with data in graph
        if tot_species < 2:
            print('- ERROR, Total number of input taxon ids must be greater than 1. For example 9606,8355... Exiting...')
            sys.exit()
        
        if len(set(t_ids) & valid_species_ids) != len(t_ids):
            print('- ERROR, One or more incompatible taxon_ids input. Exiting...')
            sys.exit()
        
    else:
        t_ids = list(org_taxon_ids.keys())
    
    tot_species = len(t_ids)
    tot_comps = (tot_species*(tot_species-1))/2
    print("- Species found for comparison {}, Total comparisons to make {}".format(tot_species, tot_comps))

    # Generate our input configurations list
    # This ensures that all filepaths exist for all comparisons before we start making any computations

    # The following double for loop performs all pairwise comparisons with out repeats (i.e. a->b but not b->a)
    # TO DO: Do we want to double this? (i.e. a->b, and b->a)
    configs = []
    for i, species_a in enumerate(t_ids[:-1]):
        a_name = ids_to_name[species_a]
        for species_b in t_ids[i+1:]:
            b_name = ids_to_name[species_b]

            # Base level config
            cpath = os.path.join(input_args.project_dir, "species_data", "common_orthologs_{}_vs_{}.tsv".format(a_name, b_name))
            apath = os.path.join(input_args.project_dir, "species_data", "{}_phenotype_to_ortholog.pkl".format(a_name))
            bpath = os.path.join(input_args.project_dir, "species_data", "{}_phenotype_to_ortholog.pkl".format(b_name))

            # Generate two configs for comparison (Species A-->B, and B-->A)
            config_a = {"species_a":a_name,
                        "species_b":b_name,
                        "species_a_phenotype_path":apath,
                        "species_b_phenotype_path":bpath,
                        "common_orthologs_path":cpath,
                        "output_directory":check_outdir}
            
            # Swap a & b values to get comparison in other direction
            config_b = {"species_a":b_name,
                        "species_b":a_name,
                        "species_a_phenotype_path":bpath,
                        "species_b_phenotype_path":apath,
                        "common_orthologs_path":cpath,
                        "output_directory":check_outdir}
            
            # Just do one config for now
            configs.append(config_a)
            #configs.append(RandomSpeciesComparison.parse_obj(config_a))
    
    print("- {} Pairwise comparisons configurations created...".format(len(configs)))
    return t_ids, configs



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
    # Basic run command (Human vs. Mouse for 100 trials across 10 cores)
    ###python phenologs_randomized_fdr_pvals.py -taxon_ids 9606,10090 -n 100 -c 10 -p path/to/top_level_project_dir/

    taxon_ids, comparison_configs = initiate_random_species_comparison_configs(args)
    for config in comparison_configs:
        print("- Computing {} random trials between {} -- {}".format(args.num_trials, 
                                                                     config["species_a"], 
                                                                     config["species_b"]))

        run_comparisons_parallel(config, n_trials=args.num_trials, num_proc=args.cpu_cores)