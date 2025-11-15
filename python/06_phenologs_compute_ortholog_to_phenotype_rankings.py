# General imports
import os
import sys
import argparse
import copy
import pandas as pd

# Custom imports
from phenologs_utils import (load_fdr_table_to_lookup,
                             pool_phenologs_data,
                             initiate_ortholog_to_phenotype_ranking_calculation_config,
                             OrthologToPhenotypeCalculations)


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
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        parser.add_argument("-fdr", help="One minus the false discovery rate.. .95 is default", required=True, type=float, default=.95)
        parser.add_argument("-kneighbs", help="k-nearest phenologs to use when combing across multiple phenologs", required=True, type=int, default=10)
        parser.add_argument("-rank_metric", help="Which metric to use for combining knearest neighbor... default is naive_bayes (nb), (hg is hyper geometric)", required=False, choices=['nb', 'hg'], default='nb')
        parser.add_argument("-rank_type", help="Which type of analyis to perform, gene ranking or protein family ranking", required=False, choices=["protein_family", "gene"], default='protein_family')
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    # Initiate a dictionary to hold filepaths and necesarry variables from our input set of arguments
    sp_config = initiate_ortholog_to_phenotype_ranking_calculation_config(args)

    # Load fdr table to lookup
    fdr_lookup = load_fdr_table_to_lookup(sp_config["fdr_path"])

    # Check if pooled phenologs file exists or not
    if os.path.isfile(sp_config["sig_phenologs_path"]):
        print("- Pooled phenologs file already exists at {}. Read into memory instead of overwiting...".format(sp_config["sig_phenologs_path"]))
        sig_phenolog_df = pd.read_csv(sp_config["sig_phenologs_path"], low_memory=False, sep='\t')
    
    else:
        # Pool significant phenologs and write to file
        sig_phenolog_df = pool_phenologs_data(sp_config["results_dir"], 
                                            fdr_lookup, 
                                            sp_config["sig_phenologs_path"], 
                                            fdr_level=sp_config["fdr"], 
                                            compress=False)
    
    # Compute gene-->phenotype rank/distance matrix 
    ##ddd = OrthologToPhenotypeCalculations.model_validate(sp_config).compute_ortholog_phenotype_distances()
    ##ddd = OrthologToPhenotypeCalculations.model_validate(sp_config).compute_ortholog_phenotype_distances()

    if args.rank_type == "protein_family":
        ddd = OrthologToPhenotypeCalculations.model_validate(sp_config).compute_ortholog_phenotype_distances()
    
    elif args.rank_type == "gene":
        ddd = OrthologToPhenotypeCalculations.model_validate(sp_config).compute_ortholog_phenotype_distances_gene_centric()
