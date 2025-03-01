# General imports
import os
import sys
import argparse
import copy
import pickle
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from collections import Counter
from typing import Dict, Optional

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
    
    print("- Total filepaths gathered {}".format(format(sum([len(v) for v in trials_by_species.values()]), ',')))
    print("- Total pairwise species comparison types found {}".format(len(trials_by_species)))
    for k, v in trials_by_species.items():
        print("- {} -- Trials found {}".format(k, format(len(v), ',')))
    
    return trials_by_species, sorted(list(trial_nums.keys()))


def collapse_cross_species_trials(trials, trials_by_species, comp_key):
    """
    Keeps track of unique p-values found across all of our trials. 
    10,000 Human phenotypes x 10,000 mouse phenotypes = 100,000,000 trial phenolog calculations
    We do 1,000 trials so we would have 100 Billion data points. Instead, we collapse our data by
    unique pvalues that are found and count the occurence of them across all trials.
    By doing this, we end up with ~941,000 key, value pairs that we are able to work with much more
    effecfively.
    """
    
    # Compute a single FDR per trial across entire species-species pairwise comparison dataset
    pvals = {}
    for trial_num, trial_path in trials_by_species[comp_key].items():
        
        # Open file, and expand pvalue list based on "hg_pval" and "occurrence" columns
        if ".gz" in trial_path: 
            df = pd.read_csv(trial_path, sep='\t', compression="gzip")
        else: 
            df = pd.read_csv(trial_path, sep='\t')
        
        # Gather pvals for this species--species comparison
        for pval, occ in zip(list(df["hg_pval"].astype(float)), list(df["occurrence"].astype(int))):
            if pval not in pvals:
                pvals.update({pval:0.})
            pvals[pval] += occ
        
        ##print(trial_num)

    print("- Total unique pvalues found {}".format(format(len(pvals), ',')))
    return pvals
    

def collapse_cross_species_phenologs(input_dir, comp_key):
    """
    Keeps track of unique p-values found across all of our trials. 
    10,000 Human phenotypes x 10,000 mouse phenotypes = 100,000,000 trial phenolog calculations
    We do 1,000 trials so we would have 100 Billion data points. Instead, we collapse our data by
    unique pvalues that are found and count the occurence of them across all trials.
    By doing this, we end up with ~941,000 key, value pairs that we are able to work with much more
    effecfively.
    """
    
    #Header column values 
    #Species A Phenotype ID  
    #Species B Phenotype ID  
    #Ortholog Count A        
    #Ortholog Count B        
    #Overlap Count   
    #Common Ortholog Count   
    #P-Value 
    #Species A Phenotype Name        
    #Species B Phenotype Name
    
    
    # Read in data to memory 
    phenologs_filepath = os.path.join(input_dir, "{}_vs_{}_all_phenologs.tsv".format(comp_key[0], comp_key[1]))
    phen_df = pd.read_csv(phenologs_filepath, sep='\t')
    
    pvals = {}
    for pval in list(phen_df["hg_pval"].astype(float)):
        if pval not in pvals:
            pvals.update({pval:0})
        pvals[pval] += 1
    
    print("- Total unique pvalues found {}".format(format(len(pvals), ',')))
    return pvals


# For visualization and FDR calculations
def collapsed_data_to_cumulitive(input_data: Dict):
    
    total = float(sum(list(input_data.values())))
    sorted_data = np.asarray(sorted([[k,v] for k,v in input_data.items()], key=lambda x: x[0], reverse=True))
    cumul = 0.
    data_frac = []
    data_counts = []
    for d in sorted_data:
        cumul += d[1]
        data_frac.append(1. - (cumul/total))
        data_counts.append(total-cumul)
    
    return sorted_data, np.asarray(data_frac), total


def get_fdr_cuttoffs(sorted_pval_counts: list):
    """
    Specific algorithm to compute the pvalues one should use to get xyz desired false discovery rate...
    Input is a specific data structure where each element within the list is [pvalue, occurence]...
    This is for a collapsed form of the data (to save memory and compute) where we can 'bin' the cumulitive
    fraction of the data to any desired level of granularity without sacrificing runtime.
    """
    
    # Algorithm to take an array where each element is the number of times that a particular value occurred
    
    # Total number of comparisons made (collapsed by pvalue occurence)
    tot = sum(np.asarray(sorted_pval_counts).T[1])
    fdrs = [i for i in np.arange(.95, 1.00001, .00001)]
    fdc = {fd:1. for fd in fdrs}

    # Alg variables
    cumul, frac = 0., 0.
    fd, fdr_count = fdrs[0], len(fdrs)
    
    # Loop through sorted [[pvalue,occurence], ...] data and compute cumulative fraction of data and the
    # corresponding pvalue at that fraction (i.e. what pvalue do we use to get a false discovery rate of xyz)
    for v in sorted_pval_counts:
        pval, occ = v
        cumul += occ
        frac = cumul/tot

        if frac >= fd:
            fdc[fd] = pval

            while True:

                if fdr_count == 0:
                    break

                fd = fdrs.pop(0)
                fdr_count -= 1
                if frac >= fd:
                    fdc[fd] = pval
                else:
                    break
    
    return fdc


def get_cumulitive_phenolog_counts(sorted_pval_counts: list, pval_cutoffs: Dict):
    """
    needs documenting
    """
    
    # Algorithm to take an array where each element is the number of times that a particular value occurred
    
    # Total number of comparisons made (collapsed by pvalue occurence)
    tot = sum(np.asarray(sorted_pval_counts).T[1])
    fdrs = sorted(list(pval_cutoffs.keys()))
    fdc = {fd:1. for fd in fdrs}

    # Alg variables
    cumul, frac = 0., 0.
    fd, fdr_count = fdrs[0], len(fdrs)
    
    # Loop through sorted [[pvalue,occurence], ...] data and compute cumulative fraction of data and the
    # corresponding pvalue at that fraction (i.e. what pvalue do we use to get a false discovery rate of xyz)
    for v in sorted_pval_counts:
        pval, occ = v
        cumul += occ
        frac = cumul/tot

        if pval <= pval_cutoffs[fd]:

            fdc[fd] = tot - cumul

            while True:

                if fdr_count == 0:
                    break

                fd = fdrs.pop(0)
                fdr_count -= 1
                if pval <= pval_cutoffs[fd]:
                    fdc[fd] = tot - cumul
                else:
                    break

    return fdc


def get_trial_fdr_data(trial_data, n_trials, results_dir):
    
    species_trial_info = {}
    for k in trial_data:
        
        # Collapse random trials, and results data to respective singular datastructures
        # NOTE - This is the step that should be parallelized!
        random_trial_pvals = collapse_cross_species_trials(trials=n_trials, trials_by_species=trial_data, comp_key=k)
        
        # Everything from here down is fine single core...
        random_bulk_info = collapsed_data_to_cumulitive(random_trial_pvals)
        random_data, random_data_frac, random_total = random_bulk_info

        # Collapse and upack real phenolog calculation data
        # This could be parallelized... but probably not worth the effort
        phenologs_pvals = collapse_cross_species_phenologs(input_dir=results_dir, comp_key=k)
        phen_bulk_info = collapsed_data_to_cumulitive(phenologs_pvals)
        phen_data, phen_data_frac, phen_total = phen_bulk_info

        # Generate false discovery rate pvalue cutoffs based on random trial nonzero intersection pvalue distribution
        fdr_cutoffs = get_fdr_cuttoffs(random_data)

        # Figure out cumulative counts of real data phenolog counts from pvalue thresholds
        fdr_cumul = get_cumulitive_phenolog_counts(phen_data, pval_cutoffs=fdr_cutoffs)

        # Add our data for second figure tallying how many phenologs above xyz fdr we get for each cross species comparison
        lab = "{}{}{} at default fdr".format(k, '\n', phen_data.T[0])
        
        species_trial_info.update({k:{"fdr_cutoffs":fdr_cutoffs, 
                                      "fdr_cumul":fdr_cumul,
                                      "label":lab,
                                      "random_data":random_bulk_info,
                                      "real_data":phen_bulk_info}})

    return species_trial_info
    

def plot_fdr_cumul_results(trial_info: Dict, fdr_default=.95, savefig1=False, savefig2=False, display_figs=False):
    
    # nrows x 2 columns formatted multiplot figure 
    nrows, ncols = int(math.ceil(len(trial_data)/2)), 2 
    fig = plt.figure()
    fig.set_size_inches(ncols*5, nrows*3) # Width x height
    
    row_ind, col_ind = 0, 0
    for k in trial_info:

        ax = plt.subplot2grid((nrows, ncols), (row_ind, col_ind))
        col_ind += 1
        if col_ind == ncols:
            col_ind = 0
            row_ind += 1

        # Random, real, fdr data
        random_data, random_data_frac, random_total = trial_info[k]["random_data"]
        phen_data, phen_data_frac, phen_total = trial_info[k]["real_data"]
        fdr_cutoffs = trial_info[k]["fdr_cutoffs"]
        
        ax.plot(random_data.T[0], random_data_frac, label="Random")
        ax.plot(phen_data.T[0], phen_data_frac, label="Real")

        # FDR data
        ymin,ymax = ax.get_ylim()
        ax.plot([fdr_cutoffs[fdr_default], fdr_cutoffs[fdr_default]], 
                [ymin, ymax], 
                linestyle='--', 
                color='red', 
                alpha=.5, 
                label="FDR {}{}p-val {}".format(round(1.-fdr_default, 4), 
                                                                     '\n', 
                                                                     "{:e}".format(fdr_cutoffs[fdr_default])))

        ax.set_xlim(min(random_data.T[0]) / 2., 1.)

        ax.invert_xaxis()
        ###ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlabel("Hypergeometric probability of observing{} gene intersection by chance".format('\n'))
        ax.set_ylabel("Fraction of phenologs")

        ax.legend()
        ax.set_title("{} vs. {}".format(k[0], k[1]))
    
    # Save and display options
    plt.tight_layout()
    if savefig1 != False:
        plt.savefig(savefig1) 
    if display_figs != False:
        plt.show()
    else:
        plt.close()
    
    ############################################################################
    ### Plot to show distribution of significant phenolog across all species ###
    fig = plt.figure()
    ax = plt.subplot2grid((1, 1), (0, 0))

    for k in trial_info:
        cdata = trial_info[k]["fdr_cumul"]
        plot_data = np.asarray(sorted([[kk,v] for kk,v in cdata.items()], key=lambda x: x[0], reverse=True)).T
        lab = " vs. ".join([k[0], k[1]])
        ax.plot(plot_data[1], plot_data[0], label=lab)

    #ax.set_ylim(fdr_default, 1.)
    ax.set_xlabel("Number of phenologs above score threshold (log10)")
    ax.set_ylabel("1 - False Discovery Rate")
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xscale("log")
    
    # Save and display options
    if savefig2 != False:
        plt.savefig(savefig2, bbox_inches='tight')
    if display_figs != False:
        plt.show()
    else:
        plt.close()


# Leveraged in phenologs_randomized_fdr_pvals
def initiate_phenologs_fdr_fpaths(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    - Input species is compared against all other species within the monarch kg
    """

    # Ensures part1 & part2 of pipeline have been completed
    check_file = os.path.join(input_args.project_dir, "species_data", "species_information_table.tsv")
    check_outdir = os.path.join(input_args.project_dir, "random_trials")

    if not os.path.isfile(check_file):
        print("- ERROR, Project species information table doesn't seem to exist. Exiting...")
        sys.exit()
    
    if not os.path.isdir(check_outdir):
        print("- ERROR, Project phenologs_results directory doesn't seem to exist. Exiting...")
        sys.exit()
    
    # Figure out which species ids / names we have and which ones are relevant
    species_df = pd.read_csv(check_file, sep='\t')
    species_df = species_df[species_df["Genes with >= 1 phenotype edge"] > 0] # Only want species whith non zero phenotype information

    # Pull out taxon_id information (twice as two separate variables)
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon label"]))}
    org_taxon_ids = copy.copy(ids_to_name)
    
    # Format, and add additional keys to make more friendly to input arguments
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon label"]))}
    ids_to_name.update({sp_id.split(":")[1]:v for sp_id,v in ids_to_name.items()}) # Add keys without NCBITaxon: prefix
    
    # Figure out which species id(s) we are tasked with comparing to one another
    valid_species_ids = set(ids_to_name.keys())
    if input_args.taxon_id == "all":
        t_ids = list(org_taxon_ids.keys())
    elif input_args.taxon_id in valid_species_ids:
        t_ids = [input_args.taxon_id]
    else:
        print('- ERROR, all or relevant taxon id must be supplied for taxon_id argument. Exiting...')
        sys.exit()
    
    
    # Now we need to filter our comparison species list by relevant / available prediction networks
    # "phenotype" or "disease" networks are available for use
    comp_species = list(org_taxon_ids.keys())
    res_paths, rand_paths, sp_names = [], [], []

    for sp_id in t_ids:
        nformatted = ids_to_name[sp_id].replace("-", "_")
        rand_name = "{}_{}_trials".format(nformatted, input_args.prediction_network)
        rand_dir = os.path.join(input_args.project_dir, "random_trials", rand_name)
        
        res_name = "{}_{}_results".format(nformatted, input_args.prediction_network)
        res_dir = os.path.join(input_args.project_dir, "phenologs_results", res_name)

        if os.path.isdir(res_dir) and os.path.isdir(rand_dir):
            res_paths.append(res_dir)
            rand_paths.append(rand_dir)
            sp_names.append(ids_to_name[sp_id])
        else:
            print("- Warning, random_trials and or phenologs_results appear to be missing for {}".format(args.taxon_id))
    
    if len(res_paths) == 0:
        print('- ERROR, No relevant data found for taxon_id {} of type {}. Exiting...'.format(args.taxon_id,
                                                                                              args.prediction_network))
        sys.exit()
    
    return res_paths, rand_paths, sp_names


# def write_fdr_tables_from_results(fdr_results, trial_table_outpath, fdr_table_outpath):
#     """
#     Writes two tables. One for the fdr averages across trials (a single rowed table minus the header),
#     and another table for fdr values found across all trials (number of rows == N trials run)
#     """
    
#     df_data = {}
#     df_fdr_data = {}
#     for k,v in fdr_results.items():
#         kk = "fdr:{}".format(k)
#         df_data.update({kk:v})
#         df_fdr_data.update({kk:[np.average(v)]})
    
#     pd.DataFrame(df_data).to_csv(trial_table_outpath, sep='\t', index=False)
#     pd.DataFrame(df_fdr_data).to_csv(fdr_table_outpath, sep='\t', index=False)
#     print("- FDR trials table file written to {}".format(trial_table_outpath))
#     print("- FDR  table file written to {}".format(fdr_table_outpath))
    

# # EXPERIMENTAL
# def read_trials_into_mem_parallel(input_dir, num_proc: int = 1):
#     """
#     Reads 
#     """

#     # Deal with - and 0 type edge cases, and instantiate all our objects before running computation
#     num_proc = max(1, num_proc)

#     # Run pipeline
#     trial_info, trials = gather_trial_data(input_dir)
#     print("- Trial data gathered...")

#     # Copy information as many times as cpu cores
#     div_trial_info = [copy.copy(trial_info) for i in range(0, num_proc)]

#     # Divide trials for parallel
#     div_trials = divide_workload(trials, num_proc=num_proc)

#     # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
#     output = mp.Queue()
#     pool = mp.Pool(processes=num_proc)
#     results = [pool.apply_async(compute_fdrs_by_trial_nums, args=(d, t)) for d, t in zip(div_trials, div_trial_info)]
#     output = [ p.get() for p in results ]

#     # Merge results from previous step into single data structure
#     # Each element in the output list is a dictionary {sival:[pval,pval,...], }
#     merged_data = {}
#     for p in output:
#         for k,v in p.items():
#             if k not in merged_data:
#                 merged_data.update({k:[]})
#             merged_data[k] += v

#     return merged_data



if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='-p argument OR -i & -o arguments are required. Data will  be pulled \
                                                      from an assumed directory / project structure. This script will gather \
                                                      all randomized trial results from speciesA-speciesB comparisons, and finds \
                                                      the pvalues corresponding to different false discovery rates. \
                                                      SpeciesA, SpeciesB, fdr, pval... for all data is written to _fdr_table.tsv')

        parser.add_argument("-p","--project_dir", help="Top most project directory", required=False, type=str, default=None)
        parser.add_argument("-c", "--cpu_cores", help="Number of cpu cores to use.", required=False, type=int, default=1)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        parser.add_argument("-display_figs", help='Option to walkthrough analysis figure by figure in separate windows', required=False, type=bool, default=False)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    fpaths, rpaths, snames = initiate_phenologs_fdr_fpaths(args)
    for res_path, rand_path, sname in zip(fpaths, rpaths, snames):
        
        # Gather random trial data and compute fdr cuttoff values from hg pvalues
        trial_data, n_trials = gather_trial_data(rand_path)
        trial_info = get_trial_fdr_data(trial_data, n_trials, res_path)
        
        # Generate our outfile paths for data and plots
        fdr_table_outpath = os.path.join(res_path, "{}_fdr_table.tsv".format(sname))
        plot_data_outpath = os.path.join(res_path, "{}_fdr_plot_data.pkl".format(sname))
        plot1_data_outpath = os.path.join(res_path, "{}_xspecies_random_vs_real.pdf".format(sname))
        plot2_data_outpath = os.path.join(res_path, "{}_xspecies_phenolog_counts.pdf".format(sname))
        
        # Generate and write fdr table
        fdr_table = np.asarray([[",".join(k),kk,v] for k in trial_info.keys() for kk,v in trial_info[k]["fdr_cutoffs"].items()]).T

        pd.DataFrame({"species_comparison":fdr_table[0],
                      "fdr":fdr_table[1], 
                      "pvalue_cutoff":fdr_table[2]}).to_csv(fdr_table_outpath, sep='\t', index=False)
        
        # Write out data used for making plots
        pickle.dump(trial_info, open(plot_data_outpath, 'wb'))
    
        # Plot our results and save figures
        plot_fdr_cumul_results(trial_info, 
                               fdr_default=.95, 
                               savefig1=plot1_data_outpath, 
                               savefig2=plot2_data_outpath,
                               display_figs=args.display_figs)
        
        print("- FDR table and figure(s) generated... Done!")