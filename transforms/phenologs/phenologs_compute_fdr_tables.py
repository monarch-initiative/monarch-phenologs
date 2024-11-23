import os
import sys
import argparse
import copy
import numpy as np
import pandas as pd
import multiprocessing as mp
from collections import Counter

# Custom imports
from phenologs_utils import divide_workload



def gather_trial_data(input_dir):
    """
    Grabs all files within the input directory where
    - "_vs_" is contained in the filename, 
    - the filename ends with ".tsv" or ".tsv.gz"
    - the last charcter of the filename before ".tsv" must pertain to a trial number (i.e. an integer)

    {(species_a, species_b):{trial_number:path/to/trial/results/file}} is the format we store all our our data
    """
    
    # Gather all filepaths to trial data
    trial_files = []
    trials_by_species = {}
    trial_nums = {}
    for fname in os.listdir(input_dir):
        
        if "_vs_" not in fname:
            continue
        
        if not fname.endswith(".gz") and not fname.endswith(".tsv"):
            continue
        
        # Ensure last piece of filename corresponds to a trial number
        fname_formatted = fname.replace(".gz", '')
        fname_formatted = fname_formatted.replace(".tsv", '')
        cols = fname_formatted.split("_")
        
        # See if we have a trial number as the last bit of our filename
        try:
            trial_num = int(cols[-1])
            trial_nums.update({trial_num:''})
        except:
            continue
        
        
        # Keep track of trial information by species-species comparison
        ddd = fname_formatted.split("_vs_")
        specA = ddd[0]
        specB = ddd[1].split("_")[0]
        key = tuple((specA, specB))

        if trials_by_species.get(key, None) == None:
            trials_by_species.update({key:{}})
        
        if trials_by_species[key].get(trial_num, None) == None:
            trials_by_species[key].update({trial_num:os.path.join(input_dir, fname)})
    
    # TO DO: Check for inconsistent number of trials between species-species comparisons.
    #        What to do if this is found? Maybe just a warning...
    print("- Total filepaths gathered {}".format(format(sum([len(v) for v in trials_by_species.values()]), ',')))
    print("- Total pairwise species comparison types found {}".format(len(trials_by_species)))
    for k, v in trials_by_species.items():
        print("- {} -- Trials found {}".format(k, format(len(v), ',')))
    
    return trials_by_species, sorted(list(trial_nums.keys()))


def compute_fdrs_by_trial_nums(trials, trials_by_species):
    """
    Predetermined false discovery rate cutoffs (i.e .05, .01, .001, .0001) are used
    to signify how far into the leading edge of the pvalue distribution we should go to grab pvalues
    to signify a line in the sand at which we call a phenolog significant or not. Creates two tables...
    The first is the pvalue cuttoff found for any given set of random trials (i.e. all trials pertaining to number 1 or 100, or 55 etc...),
    and the second table is the average of these pvalues across the total number of trials.
    """
    
    # Compute a single FDR per trial across entire species-species pairwise comparison dataset
    tot_trials = len(trials)
    trials_data = []
    sig_levels = [.05, .01, .001, .0001]
    fdr_by_sig_level = {sgl:[] for sgl in sig_levels}
    cc = 0
    for n in trials:
        pvals = []
        for k, v in trials_by_species.items():
            trial_path = v[n]

            # Open file, and expand pvalue list based on "hg_pval" and "occurrence" columns
            if ".gz" in trial_path: 
                df = pd.read_csv(trial_path, sep='\t', compression="gzip")
            else: 
                df = pd.read_csv(trial_path, sep='\t')
            
            # Gather pvals for this species--species comparison
            for pval, occ in zip(list(df["hg_pval"].astype(float)), list(df["occurrence"].astype(int))):
                pvals += ([pval] * occ)
        
        # Now that we have all pvalues compiled across all random species--species comparisons,
        # we can sort them, and come up with an "FDR" for this data set (for different significance levels)
        pvals = sorted(pvals, reverse=False) # Smallest first
        for sgl in sig_levels:
            cuttoff_pos = round((len(pvals)) * sgl)
            pval_cuttoff = pvals[cuttoff_pos]
            fdr_by_sig_level[sgl].append(pval_cuttoff)
            ### fdr_by_sig_level[sgl].append(pvals[round((len(pvals)) * sgl)]) # One liner (but way less readable)

        cc += 1
        if cc % 10 == 0:
            print("- Trial {}/{} read into memory...".format(cc, tot_trials))

    return fdr_by_sig_level


def write_fdr_tables_from_results(fdr_results, trial_table_outpath, fdr_table_outpath):
    """
    Writes two tables. One for the fdr averages across trials (a single rowed table minus the header),
    and another table for fdr values found across all trials (number of rows == N trials run)
    """
    
    df_data = {}
    df_fdr_data = {}
    for k,v in fdr_results.items():
        kk = "fdr:{}".format(k)
        df_data.update({kk:v})
        df_fdr_data.update({kk:[np.average(v)]})
    
    pd.DataFrame(df_data).to_csv(trial_table_outpath, sep='\t', index=False)
    pd.DataFrame(df_fdr_data).to_csv(fdr_table_outpath, sep='\t', index=False)
    print("- FDR trials table file written to {}".format(trial_table_outpath))
    print("- FDR  table file written to {}".format(fdr_table_outpath))
    

# EXPERIMENTAL
def read_trials_into_mem_parallel(input_dir, num_proc: int = 1):
    """
    This will distribute a list of configs to each cpu core requested (i.e. num_proc argument).
    Each config will run a single species vs species phenologs calculation and write results
    """

    # Deal with - and 0 type edge cases, and instantiate all our objects before running computation
    num_proc = max(1, num_proc)

    # Run pipeline
    trial_info, trials = gather_trial_data(input_dir)
    print("- Trial data gathered...")

    # Copy information as many times as cpu cores
    div_trial_info = [copy.copy(trial_info) for i in range(0, num_proc)]

    # Divide trials for parallel
    div_trials = divide_workload(trials, num_proc=num_proc)

    # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
    output = mp.Queue()
    pool = mp.Pool(processes=num_proc)
    results = [pool.apply_async(compute_fdrs_by_trial_nums, args=(d, t)) for d, t in zip(div_trials, div_trial_info)]
    output = [ p.get() for p in results ]

    # Merge results from previous step into single data structure
    # Each element in the output list is a dictionary {sival:[pval,pval,...], }
    merged_data = {}
    for p in output:
        for k,v in p.items():
            if k not in merged_data:
                merged_data.update({k:[]})
            merged_data[k] += v

    return merged_data



if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='-p argument OR -i & -o arguments are required. Data will either be pulled \
                                                      from an assumed directory / project structure, or the input and output \
                                                      directory paths can be supplied individually. This script will gather \
                                                      randomized trial results (on a per trial basis across all species \
                                                      and computes a false discovery rate for different significance cutoffs. \
                                                      In other words, all files pertaining to trial number 1 will have pvalues \
                                                      pooled into a single distribution, and the pvalue(s) == top5%, for example, \
                                                      would be taken as a single fdr value (we can do this for multiple cutoffs). \
                                                      So for N trials we will have N fdr values that we can derive our \
                                                      pvalue cutoff for the "final phenologs calculation". This information is \
                                                      ultimately collated into a table for ease of use downstream.')

        parser.add_argument("-p","--project_dir", help="Top most project directory", required=False, type=str, default=None)
        parser.add_argument("-c", "--cpu_cores", help="Number of cpu cores to use.", required=False, type=int, default=1)

        group1 = parser.add_argument_group()
        group1.add_argument("-i", "--input_dir", help="Path of directory containining randomized trial files", required=False, type=str, default=None)
        group1.add_argument("-o", "--output_dir", help="Path of directory to write fdr results tables to", required=False, type=str, default=None)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    # Resolve input arguments
    if args.project_dir:
        input_dir = os.path.join(args.project_dir, "random_trials")
        output_dir = os.path.join(args.project_dir, "random_trials_fdr")

    elif args.input_dir:
        input_dir = args.input_dir
        output_dir = args.output_dir
    else:
        print('- ERROR, invalid input argument combination. -p path/to/projectsdir  OR  -i path/to/trials/dir/ -o path/to/outputdir/')
        sys.exit()

    # Output filepaths
    outpath1 = os.path.join(output_dir, "fdr_trials_table.tsv")
    outpath2 = os.path.join(output_dir, "fdr_table.tsv")

    # Gather relevant filepaths, and read into memory (many files so scoop them up via multiprocess)
    trial_pvals = read_trials_into_mem_parallel(input_dir, num_proc=args.cpu_cores)

    # Write results
    write_fdr_tables_from_results(trial_pvals, outpath1, outpath2)
    print("- Done!")