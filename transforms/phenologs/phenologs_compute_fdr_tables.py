import os
import argparse
import numpy as np
import pandas as pd
from collections import Counter


def gather_trial_data(input_dir):
    
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
        specA = ddd[0].split("_")[-1]
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
    
    df_data = {}
    df_fdr_data = {}
    for k,v in fdr_results.items():
        kk = "sig:{}".format(k)
        df_data.update({kk:v})
        df_fdr_data.update({kk:[np.average(v)]})
    
    pd.DataFrame(df_data).to_csv(trial_table_outpath, sep='\t', index=False)
    pd.DataFrame(df_fdr_data).to_csv(fdr_table_outpath, sep='\t', index=False)
    print("- FDR trials table file written to {}".format(trial_table_outpath))
    print("- FDR  table file written to {}".format(fdr_table_outpath))
    
    
if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Computes randomized comparison trials between species a vs species b')
        parser.add_argument("-i", "--input_dir", help="Path of directory containining randomized trial files", required=True, type=str)
        parser.add_argument("-o", "--output_dir", help="Path of directory to write fdr results tables to", required=True, type=str)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    outpath1 = os.path.join(args.output_dir, "fdr_table.tsv")
    outpath2 = os.path.join(args.output_dir, "fdr_trials_table.tsv")

    trial_info, trials = gather_trial_data(args.input_dir)
    trial_pvals = compute_fdrs_by_trial_nums(trials=trials, trials_by_species=trial_info)
    write_fdr_tables_from_results(trial_pvals, outpath1, outpath2)


