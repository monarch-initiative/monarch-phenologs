# General imports
import os
import sys
import argparse
import copy
import pandas as pd

# Custom imports
from phenologs_utils import (load_fdr_table_to_lookup,
                             pool_phenologs_data,
                             OrthologToPhenotypeCalculations)


def initiate_ortholog_to_phenotype_ranking_calculation_config(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    - Input species is compared against all other species within the monarch kg
    """

    # Ensures part1 & part2 of pipeline have been completed
    species_data_dir = os.path.join(input_args.project_dir, "species_data")
    check_file = os.path.join(species_data_dir, "species_information_table.tsv")
    check_outdir = os.path.join(input_args.project_dir, "phenologs_results")

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
    if input_args.taxon_id in valid_species_ids:
        sp_id = input_args.taxon_id
    else:
        print('- ERROR, relevant taxon id must be supplied for taxon_id argument. Exiting...')
        sys.exit()
    
    # Now we need to filter our comparison species list by relevant / available prediction networks
    # "phenotype" or "disease" networks are available for use
        
    # Format filenames and paths for fdr table
    sp_name = ids_to_name[sp_id]
    nformatted = ids_to_name[sp_id].replace("-", "_")
    res_name = "{}_{}_results".format(nformatted, input_args.prediction_network)
    res_dir = os.path.join(input_args.project_dir, "phenologs_results", res_name)
    fdr_path = os.path.join(res_dir, "{}_fdr_table.tsv".format(sp_name))
    
    
    # Check if necessary fdr table exists and gather remaining filepaths necessary
    if os.path.isdir(res_dir) and os.path.isfile(fdr_path):

        # Our species specific table filepaths from initial steps of pipeline
        p2o_path = os.path.join(species_data_dir, "{}_{}_to_ortholog.pkl".format(sp_name, input_args.prediction_network))
        g2o_path = os.path.join(species_data_dir, "{}_gene_to_ortholog.tsv".format(sp_name))
        g2p_path = os.path.join(species_data_dir, "{}_gene_to_{}.tsv".format(sp_name, input_args.prediction_network))
        
        # Pooled phenolog file name (significance cutoff included in name)
        sig_phenologs_outname = "{}_pooled_phenologs_fdr{}.tsv".format(sp_name, input_args.fdr)
        sig_phenologs_outpath = os.path.join(res_dir, sig_phenologs_outname)
        
        
        # Generate dictionary to hold our filepaths for easy processing downstream
        sp_file_info = {"project_dir":input_args.project_dir,
                        "results_dir":res_dir,
                        "taxon_id":input_args.taxon_id,
                        
                        "prediction_network":input_args.prediction_network,
                        "fdr_path":fdr_path,
                        "fdr":input_args.fdr,
                        "kneighbs":input_args.kneighbs,
                        "rank_metric":input_args.rank_metric,
                        "sig_phenologs_path":sig_phenologs_outpath,
                        
                        "phen_to_orth_path":p2o_path,
                        "gene_to_orth_path":g2o_path,
                        "gene_to_phen_path":g2p_path,
                        "species_name":sp_name} # Note, phen_path can be disease or phenoptype (but phen is the general programmtic name)
        
    return sp_file_info


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
        parser.add_argument("-rank_metric", help="Which metric to use for combining knearest neighbor... default is naive_bayes", required=False, default='naive_bayes')
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
    ####if os.path.isfile(sp_config["sig_phenologs_path"]):
    ####    sig_phenolog_df = pd.read_csv(sp_config["sig_phenologs_path"])

    # Pool significant phenologs and write to file
    sig_phenolog_df = pool_phenologs_data(sp_config["results_dir"], 
                                          fdr_lookup, 
                                          sp_config["sig_phenologs_path"], 
                                          fdr_level=sp_config["fdr"], 
                                          compress=False)

    # Compute gene-->phenotype rank/distance matrix 
    ddd = OrthologToPhenotypeCalculations.model_validate(sp_config).compute_ortholog_phenotype_distances()