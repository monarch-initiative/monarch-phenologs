# General imports
import os
import sys
import shutil
import argparse
import copy
import pandas as pd
import multiprocessing as mp
from pathlib import Path

# Custom imports
from phenologs_utils import (OrthologToPhenotypeCalculations,
                             initiate_ortholog_to_phenotype_ranking_calculation_config)


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='OLD DOCUMENTATION NEEDS FIXING \
                                                     -p argument OR -i & -o arguments are required. Data will either be pulled \
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
        parser.add_argument("-c", "--cpu_cores", help="Number of cpu cores to use. HPC or long runtime is likely required", required=False, type=int, default=1)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        parser.add_argument("-fdr", help="One minus the false discovery rate.. .95 is default", required=True, type=float, default=.95)
        parser.add_argument("-kneighbs", help="k-nearest phenologs to use when combing across multiple phenologs", required=True, type=int, default=10)
        parser.add_argument("-rank_metric", help="Which metric to use for combining knearest neighbor... default is naive_bayes", required=False, default='naive_bayes')
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    num_proc = max(1, args.cpu_cores)

    # Initiate a dictionary to hold filepaths and necesarry variables from our input set of arguments
    sp_config = initiate_ortholog_to_phenotype_ranking_calculation_config(args)

    # Find the leave one out xvalidation directory, and grab relevant filepaths 
    species_name = sp_config["species_name"]
    xvalid_dir = os.path.join(args.project_dir, 
                              "leave_one_out_xvalidation", 
                              "{}_{}_results".format(species_name, args.prediction_network))
    
    # Make sure our relevant validation phenologs results directory has been made
    if not os.path.isdir(xvalid_dir):
        print("- Error, No leave one out phenologs reults found at {}... Exiting".format(xvalid_dir))
        sys.exit()

    # Grab our relevant results filepaths
    rel_suffx1 = "_pooled_phenologs_fdr{}.tsv.gz".format(args.fdr)
    rel_suffx2 = "_pooled_phenologs_fdr{}.tsv".format(args.fdr)
    res_paths = [os.path.join(xvalid_dir, fname) for fname in os.listdir(xvalid_dir) if fname.endswith(rel_suffx1) or 
                                                                                        fname.endswith(rel_suffx2)]
    # Make sure we have files to process
    if len(res_paths) == 0:
        print("- Error, No pooled_phenologs reults found at {}... Exiting".format(xvalid_dir))
        sys.exit()

    # Copy and batch base config
    div_process_objs = []
    for i in range(0, len(res_paths)):
        new_config = copy.copy(sp_config)

        # Need to alter the results file path that it has by default so we can process it parallel
        new_config["sig_phenologs_path"] = res_paths[i]
        div_process_objs.append(OrthologToPhenotypeCalculations.model_validate(new_config))
    
    # Compute gene-->phenotype rank/distance matrix 
    ##ddd = OrthologToPhenotypeCalculations.model_validate(sp_config).compute_ortholog_phenotype_distances(sp_config["sig_phenologs_path"])
    print("- Processing {} pooled results files...".format(format(len(div_process_objs), ',')))
    output = mp.Queue()
    pool = mp.Pool(processes=num_proc)
    results = [pool.apply_async(orth_obj.compute_ortholog_phenotype_distances, args=()) for orth_obj in div_process_objs]
    output = [ p.get() for p in results ]
    pool.close()
    print("- Done!")