# General imports
import os
import sys
import pickle
import random
import argparse
import copy
import time ### Can remove in the future
import numpy as np
import pandas as pd
import multiprocessing as mp
from scipy.stats import hypergeom
from collections import Counter
from IPython.display import display

from pydantic import BaseModel
from typing import Union, Literal, Dict, List, Any, Optional
import resource

# Custom imports
from phenologs_utils import species_dict


def divide_workload(data_list, num_proc: int=1) -> list:
    """
    Meant to divide up the elements in data_list into num_proc equal portions
    by iteratively adding each element to a basket never repeating the same basket until all baskets have an equal amount
    If num_proc == 1 then the original input list will be returned nested in a top layer list i.e. [data_list]
    """

    # Deal with our edge case at the very begginning which then is used as input into the second potential edge case
    ndata_elements = len(data_list)
    if ndata_elements < num_proc:
        num_proc = ndata_elements

    # Edge case
    if num_proc <= 1:
        return [data_list]
    else:
        baskets = [[] for i in range(0, num_proc)]
        index_count = 0
        for d in data_list:
            baskets[index_count].append(d)
            if index_count == (num_proc-1):
                index_count = 0
            else:
                index_count += 1

        #print("- Workload divided into {} portions with each portion recieving {} elements respectively...".format(num_proc, [format(len(b), ',') for b in baskets]))
        return baskets


def bulk_compute_hyper_geom(params):
    return {pr:float(hypergeom.pmf(*pr)) for pr in params}


class RandomSpeciesComparison(BaseModel):
        
    species_a: str
    species_b: str
    species_a_phenotype_path: str
    species_b_phenotype_path: str
    common_orthologs_path: str

    base_a_p2g: Optional[dict] = None
    base_b_p2g: Optional[dict] = None
    base_pool: Optional[list] = None
    common_ortholog_count: Optional[int] = None
    output_directory: Optional[str] = None


    def load_species_ortholog_information(self):
        """
        Initiates necessary attributes for downstream computations.
        - Number of common orthologs
        - The common ortholog pool to sample from for randomized data
        - Base phenotype-->ortholog dictionaries for each species
        """
        
        # Load species dict and common orthologs and set remaining filepaths
        base_a_p2g = pickle.load(open(self.species_a_phenotype_path, 'rb'))
        base_b_p2g = pickle.load(open(self.species_b_phenotype_path, 'rb'))
        corth_df = pd.read_csv(self.common_orthologs_path, sep='\t', header=0, low_memory=False)
        common_orthologs = list(corth_df["ortholog_id"])
        common_ortholog_count = len(common_orthologs)
        
        # Set our reusable attributes
        self.base_a_p2g = base_a_p2g
        self.base_b_p2g = base_b_p2g
        self.base_pool = common_orthologs
        self.common_ortholog_count = common_ortholog_count

        print("- Species {} vs. {} comparison data loaded... {} common orthologs found".format(self.species_a,
                                                                                               self.species_b,
                                                                                               format(self.common_ortholog_count, ',')))
    

    def run_randomized_comparison_pipeline(self, seed: Optional[int] = None):
        """
        Generates a randomized phenotype-->ortholog dataset for each species.
        Pairwise comparisons between inter species phenotypes are made
        by computing the number of orthologs that overlap between them. Note that
        ths does NOT compute pvalues, but rather gathers the set(s) of parameters used
        to calculate pvalues and full dataset of non zero overlap count parameters.
        """
        
        # Optional random seed
        if seed:
            np.random.seed(seed)
        
        # Load initial data
        self.load_species_ortholog_information()
        #print(self.species_a, self.species_a_total_genes, max([len(v) for v in self.base_a_p2g.values()]))
        #print(self.species_b, self.species_b_total_genes, max([len(v) for v in self.base_b_p2g.values()]))
        


        # NEW - We can subset out only the common orthologs for each species
        common_orths = set(self.base_pool)
        a_sub, b_sub = {}, {}
        for k,v in self.base_a_p2g.items():
            orth_list = list(set([vv for vv in v if vv in common_orths]))
            if len(orth_list) > 0:
                a_sub.update({k:orth_list})
        
        for k,v in self.base_b_p2g.items():
            orth_list = list(set([vv for vv in v if vv in common_orths]))
            if len(orth_list) > 0:
                b_sub.update({k:orth_list})

        #print("- {}, {}".format(len(a_sub), len(b_sub)))
        #print("- {}".format(self.common_ortholog_count))
        #print("- {}, {}".format(max([len(v) for v in a_sub.values()]), max([len(v) for v in b_sub.values()])))

        # Generate randomized dataset
        randomized_a = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in a_sub.items()}
        randomized_b = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in b_sub.items()}

        # OLD WAY
        # Generate randomized dataset
        ##randomized_a = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in self.base_a_p2g.items()}
        ##randomized_b = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in self.base_b_p2g.items()}
        #print("- Randomized data generation complete")
        
        # Convert to compact list data structure (instead of big numpy array)
        b_matrix = [v for v in randomized_b.values()]
        b_matrix_lengths = {i:len(v) for i, v in enumerate(b_matrix)}
        hg_a_count_params = {k:len(v) for k,v in randomized_a.items()}
        hg_b_count_params = {k:len(v) for k,v in randomized_b.items()}
        print("- Orthologs converted to integers and compact list...")

        # Map each unique value (ortholog id) in our data to the rows it belongs to (pre processing step)
        orth_to_coords = {}
        for i, row in enumerate(b_matrix):
            for v in row:
                if v not in orth_to_coords:
                    orth_to_coords.update({v:[]})
                orth_to_coords[v].append(i)
        
        #print("- Total rows {}".format(format(len(orth_to_coords), ',')))
        #print("- Total row associations {}".format(format(sum([len(v) for v in orth_to_coords.values()]), ',')))
        
        # Precompute b matrix length for initiation of zeros array within for loop
        b_phenotype_count = len(b_matrix)
        processed = 0
        t_comps = 0
        start_time = time.time()
        
        # For keeping track of hour hg params
        hg_a_counts = []
        hg_b_counts = []
        hg_overlap_counts = []
        hg_world_counts = []
        
        # Loop through all phenotype-ortholog associatoins for species a, and compute overlap with species b 
        for a_phen, a_orthologs in randomized_a.items():
            
            # Initialize datastructure to count the number of overlaps found for each of species b phen-orthologs
            counts = np.zeros(b_phenotype_count)
            a_ortho_count = hg_a_count_params[a_phen]
            
            # Compute commonality between each ortholog and species b orthologs
            for orth in a_orthologs:
                if orth in orth_to_coords:
                    counts[orth_to_coords[orth]] += 1
                t_comps += 1
            
            # This gives back an oddly shaped data array that we can flatten to get all non-zero row indices
            inds = np.argwhere(counts > 0).flatten()
            ind_count = len(inds)
            

            # TO DO: PROPER HG PARAMETERS NEED TO BE FIGURED OUT. The general framework is here, but the exact params need fixing
            # DEFAULT / PREVIOUS WAY... 
            # This will break when the number of common orthologs is < than the total number of genes assocaited with a phenotype
            # Update our hg parameter lists
            hg_overlap_counts += list(counts[inds])
            hg_b_counts += [b_matrix_lengths[ind] for ind in inds]
            hg_a_counts += [a_ortho_count] * ind_count
            hg_world_counts += [self.common_ortholog_count] * ind_count

            # Alters hg world counts parameter so that the "world size" is derived from species a total genes
            # (the actual pool of genes we are pulling from) rather than the common set of orthologs 
            #hg_overlap_counts += list(counts[inds])
            #hg_b_counts += [b_matrix_lengths[ind] for ind in inds]
            #hg_a_counts += [a_ortho_count] * ind_count
            #hg_world_counts += [self.species_a_total_genes] * ind_count
                        
            processed += 1
            if processed % 1000 == 0:
                stop_time = time.time()
                #print("- Processed {}/10,000 -- 10 samples processed in {}".format(processed, stop_time-start_time))
                #print("- {}".format(counts.shape))
                start_time = time.time()
                
        
        
        
        
        # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
        # c = number of common orthologs between phenotypes (ortholog matches)
        # N = total number of orthologs shared between species
        # m = number of orthologs in species B phenotype
        # n = number of orthologs in species A phenotype
        
        # Combine our data into a list of hg params
        tt = np.asarray([hg_overlap_counts, hg_world_counts, hg_b_counts, hg_a_counts]).T.astype(int).tolist()
        #print("- Total non zero hg params computed {}".format(format(len(tt), ',')))
        
        # Create unique set of tuples
        hg_formatted = list(map(tuple, tt))
        hg_params = set(hg_formatted)
        #print("- Unique hg params computed {}".format(format(len(hg_params), ',')))


        # Deal with occurence where the total number of orthologs associated with a given gene is more than the total common
        # orthologs between the two species



        
        # TO DO: Reduce output file format here? (i.e w)
        # Output file compression here (instead of writing every single non-zero param, we can collapse)
        # This should save ~10 fold number of lines written for each file
        tt = Counter(hg_formatted)

        return hg_params, tt
        

    def run_randomized_comparison_trials(self, n_trials: List[int]):
        
        # Keep track of unique sets of hyper geometric tests to compute so we don't recalculate
        hg_params = set()
        data = {}
        for i, n_trial in enumerate(n_trials):
            # Generate randomized data
            trial_hg_params, trial_data = self.run_randomized_comparison_pipeline()
            
            # Update our hire level output datastructures
            hg_params = hg_params | trial_hg_params
            data.update({n_trial:trial_data})
            print("- Total hg_params computed {} for {}/{} trials".format(format(len(hg_params), ','), 
                                                                          i+1, 
                                                                          format(len(n_trials), ',')))
        
        print("- {} vs. {} pairwise phenotype-->ortholog networks overlaps computed for {} trials".format(self.species_a,
                                                                                                          self.species_b,
                                                                                                          len(n_trials)))
        
        # Now compute our hyper geometric tests in bulk
        hg_pvals = bulk_compute_hyper_geom(hg_params)
        print("- Hypergeometric tests computed {}...".format(format(len(hg_pvals), ',')))
        
        # Map to pvalues, and write data files
        for trial_num, trial_data in data.items():
            
            # Define our outfile path and map params to pvalues 
            outfile_path = os.path.join(self.output_directory, 
                                        "{}_vs_{}_{}.tsv.gz".format(self.species_a, self.species_b, trial_num))
                                        
            # Map data back to computed pvalues, and format our trial data for easy writing to file
            pvals = [hg_pvals[tuple(k)] for k in trial_data]
            #param_data = np.asarray(trial_data).T
            
            # Reduced file size version
            param_data = np.asarray(list(trial_data.keys())).T
            occurrence = list(trial_data.values())
            
            # Write data using pandas
            pd.DataFrame({"a_ortholog_count":param_data[0],
                          "b_ortholog_count":param_data[1],
                          "overlap_count":param_data[2],
                          "common_ortholog_count":param_data[3],
                          "hg_pval":pvals, ## Full file size version).to_csv(outfile_path, sep='\t', index=False)

                        # Compact / reduced file size version (10 fold line count reduction for human mouse comparison)
                        "occurrence":occurrence}).to_csv(outfile_path, sep='\t', index=False, compression='gzip')

            print("- Data written to file for trial {}...".format(trial_num))
        
        ##return hg_pvals, data
        return




    # EXPERIMENTAL (Need a separate function to write data for run_comparisons_parallel_v2)
    def write_trial_data(self, trial_list, hg_pvals):
        
        for data in trial_list:

            # Map to pvalues, and write data files
            for trial_num, trial_data in data.items():
                
                # Define our outfile path and map params to pvalues 
                outfile_path = os.path.join(self.output_directory, 
                                            "{}_vs_{}_{}.tsv.gz".format(self.species_a, self.species_b, trial_num))
                                            
                # Map data back to computed pvalues, and format our trial data for easy writing to file
                pvals = [hg_pvals[tuple(k)] for k in trial_data]
                #param_data = np.asarray(trial_data).T
                
                # Reduced file size version
                param_data = np.asarray(list(trial_data.keys())).T
                occurrence = list(trial_data.values())
                
                # Write data using pandas
                pd.DataFrame({"a_ortholog_count":param_data[0],
                            "b_ortholog_count":param_data[1],
                            "overlap_count":param_data[2],
                            "common_ortholog_count":param_data[3],
                            "hg_pval":pvals, ## Full file size version).to_csv(outfile_path, sep='\t', index=False)

                            # Compact / reduced file size version (10 fold line count reduction for human mouse comparison)
                            "occurrence":occurrence}).to_csv(outfile_path, sep='\t', index=False, compression='gzip')

                print("- Data written to file for trial {}...".format(trial_num))
        
        ##return hg_pvals, data
        return


    # EXPERIMENTAL This version is more piece mail and will not compute hg pvals or write data
    def run_randomized_comparison_trials_v2(self, n_trials: List[int]):
        
        # Keep track of unique sets of hyper geometric tests to compute so we don't recalculate
        hg_params = set()
        data = {}
        for i, n_trial in enumerate(n_trials):
            # Generate randomized data
            trial_hg_params, trial_data = self.run_randomized_comparison_pipeline()
            
            # Update our hire level output datastructures
            hg_params = hg_params | trial_hg_params
            data.update({n_trial:trial_data})
            print("- Total hg_params computed {} for {}/{} trials".format(format(len(hg_params), ','), 
                                                                          i+1, 
                                                                          format(len(n_trials), ',')))
        
        print("- {} vs. {} pairwise phenotype-->ortholog networks overlaps computed for {} trials".format(self.species_a,
                                                                                                          self.species_b,
                                                                                                          len(n_trials)))
        
        return data, hg_params
        


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
    # Integers will automatically be cast to list of integers according to trial number
    if type(n_trials) == type(1):
        n_trials = [i for i in range(0, n_trials)]

    # Evenly (as possible) divide our data into baskets within a basket (list[list,list,list,...])
    if len(n_trials) > 1:
        div_trials = divide_workload(n_trials, num_proc=num_proc)
    else:
        div_trials = [n_trials]
    
    # Instantiate all our objects before running computation
    run_objs = [RandomSpeciesComparison.parse_obj(config)
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
    output = [ p.get() for p in results ]
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
    check_dir = os.path.join(input_args.project_dir, "random_trials")
    check_file = os.path.join(input_args.project_dir, "species_data", "species_information_table.tsv")
    check_outdir = os.path.join(input_args.project_dir, "random_trials_results")

    if not os.path.isdir(check_dir):
        print("- ERROR, Project directory {}/random_trials does not seem to exist. Exiting...".format(input_args.project_dir))
        sys.exit()
    
    if not os.path.isfile(check_file):
        print("- ERROR, Project species information table doesn't seem to exist. Exiting...")
        sys.exit()
    
    if not os.path.isdir(check_outdir):
        print("- ERROR, Project random_trials_results directory doesn't seem to exist. Exiting...")
        sys.exit()
    
    # Figure out which species ids / names we have and which ones are relevant
    species_df = pd.read_csv(check_file, sep='\t')
    species_df = species_df[species_df['Total Phenotype Edges'] > 0] # Only want species whith non zero phenotype information
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon Label"]))}
    ids_to_name.update({sp_id.split(":")[1]:v for sp_id,v in ids_to_name.items()}) # Add keys without NCBITaxon: prefix
    display(species_df)

    # Simple table mapping species name to the total number of genes (with at least one phenotype term associated with them)
    species_g_count = {t:int(v) for t,v in zip(list(species_df["Taxon Label"]), list(species_df["Genes >= 1 Phenotype"]))}

    # Figure out which species ids we are tasked with comparing to one another
    valid_species_ids = set(ids_to_name.keys())
    if not args.all:
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
        t_ids = list(ids_to_name.keys())
    
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
        parser.add_argument("-c", "--cpu_cores", help="Number of cores to use.", required=False, type=int, default=1)
        parser.add_argument("-taxon_ids", help="Comma separated list of taxon_ids to use", required=False, type=str)
        parser.add_argument("-all", help="Compare all available species to one another", required=False, action='store_true')
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###
    # Basic run command (Human vs. Mouse for 100 trials across 10 cores)
    ###python phenologs_randomized_fdr_pvals.py -taxon_ids 9606,10090 -n 100 -c 10 -p path/to/top_level_project_dir/

    taxon_ids, comparison_configs = initiate_pairwise_comparison_configs(args)
    for config in comparison_configs:
        print("- Computing {} random trials between {} -- {}".format(args.num_trials, 
                                                                     config["species_a"], 
                                                                     config["species_b"]))

        run_comparisons_parallel(config, n_trials=args.num_trials, num_proc=args.cpu_cores)