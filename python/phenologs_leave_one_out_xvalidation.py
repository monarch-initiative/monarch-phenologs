# General imports
import os
import sys
import pickle
import shutil
import argparse
import copy
import pandas as pd
import multiprocessing as mp
from pathlib import Path

# Custom imports
from phenologs_utils import (divide_workload,
                             load_fdr_table_to_lookup,
                             pool_phenologs_data,
                             initiate_phenologs_species_comparison_configs,
                             PhenologsSpeciesComparison)


def initiate_species_specific_filepaths(input_args):
    """
    Pulls filepaths from the project_dir/species_data for input species
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
                        "sig_phenologs_path":sig_phenologs_outpath,
                        "phen_to_orth_path":p2o_path,
                        "gene_to_orth_path":g2o_path,
                        "gene_to_phen_path":g2p_path,
                        "species_name":sp_name} # Note, phen_path can be disease or phenoptype (but phen is the general programmtic name)
        
    return sp_file_info


    
        
        # Generate dictionary to hold our filepaths for easy processing downstream


def batch_compute_leave_one_out(batch_configs, fdr_lookup_table):
    
    # Load fdr data into memory

    # Load basic config for last portion of pipeline (ortholog_to_phenotype) so we only have to do it once

    # Each config_set is a leave xyz out comparison set that we need to make (multiple species comparisons)
    for config_set in batch_configs:
        
        # Pull out variables from config and create new ones
        species_name = config_set[0]["species_a"]
        outdir = config_set[0]["output_directory"]
        outdir_name = str(Path(outdir).name)
        parent_outdir = str(Path(outdir).parent)
        fdr_level = config_set[0]["leave_out_validate_fdr"]
        
        # Make necessary output directory
        if not os.path.isdir(outdir):
            os.makedirs(outdir, exist_ok=True)
 
        # Each config is a speciesA to speciesB comparison that we need to calculate
        for config in config_set:
            PhenologsSpeciesComparison.model_validate(config).compute_cross_species_phenologs()
        
        # Pool together _all_phenologs.tsv and filter via fdr pvalue (and write results)
        # First, format pooled output filename / path
        sig_outname = "{}_{}_pooled_phenologs_fdr{}.tsv.gz".format(species_name, outdir_name, fdr_level)
        sig_outpath = os.path.join(parent_outdir, sig_outname)
        
        # Next, gather up the filepaths of the files we just wrote
        # And pool data, filter by fdr values from table, and write data
        sig_phenologs_df = pool_phenologs_data(outdir, fdr_lookup_table, sig_outpath, fdr_level, compress=True)

        # Note, we don't compute the rest of the pipeline here, so we can test different k-nearsest phenologs
        # in parallel as a separate script (so this is run first once, and then for the different k nearest neighbs)
        # we read in these results and compute

        # Delete all data within leave one out directory we just computed so we don't
        # create massive amount of data (~7,500 datasets we would end up making for human)
        shutil.rmtree(outdir)
    
    return
        

if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Leave one gene out cross validation for phenologs calculations')
        parser.add_argument("-p","--project_dir", help="Top most project directory", required=False, type=str, default=None)
        parser.add_argument("-c", "--cpu_cores", help="Number of cpu cores to use. HPC or long runtime is likely required", required=False, type=int, default=1)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        parser.add_argument("-fdr", help="One minus the false discovery rate.. .95 is default", required=True, type=float, default=.95)
        return parser.parse_args()

    args = parse_input_command()
    ############################
    
    ###############
    ### PROGRAM ###

    num_proc = max(1, args.cpu_cores)

    # Grabs relevant fdr_path (overkill at the moment), and species name 
    sp_config = initiate_species_specific_filepaths(args)

    # Load fdr data and make copies for multiprocess
    fdr_lookup = load_fdr_table_to_lookup(sp_config["fdr_path"])
    div_fdr_lookups = [copy.copy(fdr_lookup) for i in range(0, num_proc)]

    # Inititate set of base level configurations to serve as base and                    
    t_ids, base_configs = initiate_phenologs_species_comparison_configs(args)

    # Gather our "global" set of orthologs shared by speciesA and at least one other species
    global_orths = list(set([panther_id for c in base_configs for panther_id in list(pd.read_csv(c["common_orthologs_path"], sep='\t')["ortholog_id"])]))
    print("- {} global common ortholog set found".format(format(len(global_orths), ',')))

    # Load our disease/phenotype_to_ortholog data
    p2o = pickle.load(open(sp_config["phen_to_orth_path"], 'rb'))
    for k in p2o:
        p2o[k] = set(p2o[k]) # Collapse from list to set so we can lookup easier

    # Now trim back the set of orthologs to only those that have at least one connection to a disease/phenotype
    # We can only make statements about phenotypes that have more that one neighboring gene/ortholog
    # Therefore, if an ortholog's only connection is to one of these types of phenotypes it will also be removed.
    # This is ideal, because it means we do not have to perform a leave one out validation on these orthologs.
    relative_orthologs = {oid:'' for orth_set in p2o.values() for oid in orth_set if (oid in global_orths) and (len(orth_set) >= 2)}
    print("- {} orthologs found with >= 1 phenotype assocation and found within global common ortholog set".format(format(len(relative_orthologs), ',')))

    # Make base level xvalidate "results" directory 
    species_name = base_configs[0]["species_a"]
    xvalid_dir = os.path.join(args.project_dir, 
                              "leave_one_out_xvalidation", 
                              "{}_{}_results".format(species_name, args.prediction_network))

    if not os.path.isdir(xvalid_dir):
        os.makedirs(xvalid_dir, exist_ok=True)


    # Create set of configs for each ortholog we need to remove from our dataset by 
    # creating a copy, and then altering the copy's relevant key,value pairs
    xvalidate_config_sets = []
    for orth_id in relative_orthologs:
        
        # This directory will ultimitaly be made, then deleted by the same process
        process_dir = os.path.join(xvalid_dir, orth_id)
        
        packaged_configs = []
        for bconf in base_configs:
            bconf_copy = copy.copy(bconf)
            bconf_copy["output_directory"] = process_dir
            bconf_copy["leave_out_validate_set"] = set([orth_id])
            bconf_copy["leave_out_validate_fdr"] = float(args.fdr)
            packaged_configs.append(bconf_copy)
        
        # Add to our top level config holder
        xvalidate_config_sets.append(packaged_configs)
        

    # Divy up our xvalidate datasets to calculate
    ##div_configs = divide_workload(xvalidate_config_sets[0:20], num_proc=num_proc) # Limiting to 20 for testing
    div_configs = divide_workload(xvalidate_config_sets, num_proc=num_proc)

    # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
    output = mp.Queue()
    pool = mp.Pool(processes=num_proc)
    results = [pool.apply_async(batch_compute_leave_one_out, args=(dvc, fdcl)) for dvc,fdcl in zip(div_configs, div_fdr_lookups)]
    output = [ p.get() for p in results ]
    pool.close()
    print("- Done!")