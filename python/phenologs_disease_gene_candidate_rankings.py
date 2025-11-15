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
        parser = argparse.ArgumentParser(description="Script for writing gene-disease predictions across all \
                                         available phenologs data using a defined set of parameters presumably used to \
                                         process all previous data. Will write all .tsv files for each mondo \
                                         term found with phenolog data. Note, mendelian diseases are also included \
                                         whereby the actual disease causing gene is not actually recovered in the analysis. \
                                         This is because phenologs are initially computed across protein families, where \
                                         the species in comparison does not necessarily have to have the direct human ortholog \
                                         associated with the phenotype.")

        parser.add_argument("-p","--project_dir", help="Top most project directory", required=False, type=str, default=None)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        parser.add_argument("-fdr", help="One minus the false discovery rate.. .95 is default", required=True, type=float, default=.95)
        parser.add_argument("-kneighbs", help="k-nearest phenologs to use when combing across multiple phenologs", required=True, type=int, default=10)
        parser.add_argument("-rank_metric", help="Which metric to use for combining knearest neighbor... default is naive_bayes (nb), (hg is hyper geometric)", required=False, choices=['nb', 'hg'], default='nb')
        parser.add_argument("-rank_type", help="Which type of analyis to perform, gene ranking is only supported for this analysis", required=False, choices=["gene"], default='gene')
        return parser.parse_args()

    args = parse_input_command()
    ############################

    ###############
    ### PROGRAM ###

    # Initiate a dictionary to hold filepaths and necesarry variables from our input set of arguments
    sp_config = initiate_ortholog_to_phenotype_ranking_calculation_config(args)

    # Format our base output directory here
    base_outdir = os.path.join(sp_config["project_dir"], 
                              "phenologs_results", 
                              "{}_{}_results".format(sp_config["species_name"].replace("-", "_"),
                                                     args.prediction_network))

     # Format our input filenames here
    sig_phenologs_resname = "{}_pooled_phenologs_fdr{}.tsv".format(sp_config["species_name"],args.fdr)
    outdir_name = "{}_{}_to_{}_candidates_fdr{}_{}kneighbs_{}".format(sp_config["species_name"].replace("-", "_"),
                                                                      args.rank_type,
                                                                      args.prediction_network,
                                                                      args.fdr, 
                                                                      args.kneighbs, 
                                                                      args.rank_metric)
    
    g2p_resname = "{}_{}_to_{}_{}kneighbs_{}_{}.pkl".format(sp_config["species_name"],
                                                     args.rank_type,
                                                     args.prediction_network,
                                                     args.kneighbs,
                                                     args.rank_metric,
                                                     args.fdr)
    
    # Combine into paths
    outdir = os.path.join(base_outdir, outdir_name)
    sig_phens_path = os.path.join(base_outdir, sig_phenologs_resname)
    g2p_path =  os.path.join(base_outdir, g2p_resname)

    # Ensure our input / output data exists or is able to be created
    if not os.path.isfile(g2p_path):
        raise FileNotFoundError("- precomputed gene to phenotype distance file not found {}...".format(g2p_path))

    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)

    # Write our disease gene rankings
    ddd = OrthologToPhenotypeCalculations.model_validate(sp_config).write_disease_gene_rankings(outdir, g2p_path, sig_phens_path=sig_phens_path)
