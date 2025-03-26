# General imports
import os
import sys
import pickle
import random
import copy
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.stats import pearsonr
from collections import Counter
from IPython.display import display
from pathlib import Path
from pydantic import BaseModel
from typing import Dict, List, Optional


###################################################################################################
### Basic species information building blocks, and file path validation for downstream analysis ###
class Species:
    def __init__(self, taxon_id, species_name, project_dir):
        self.taxon_id = taxon_id
        self.species_name = species_name
        self.project_dir = project_dir
        self.file_exts = ["disease_to_ortholog.pkl", 
                          "gene_to_disease.tsv", 
                          "gene_to_ortholog.tsv", 
                          "gene_to_phenotype.tsv", 
                          "phenotype_to_ortholog.pkl"]

        self.disease_to_ortholog_path = None
        self.gene_to_disease_path = None
        self.gene_to_ortholog_path = None
        self.gene_to_phenotype_path = None
        self.phenotype_to_ortholog_path = None
        self.common_orthologs_paths = []
    
    def initiate_filepaths(self):
        
        # Pull in species specific file paths
        for fext in self.file_exts:
            fname = "{}_{}".format(self.species_name, fext)
            fpath = os.path.join(self.project_dir, "species_data", fname)
            if os.path.isfile(fpath):
                setattr(self, "{}_path".format(fext.split(".")[0]), fpath)
            elif "disease" not in fext:
                print("- ERROR, Species {} is missing file: {}".format(self.species_name, fext))
                sys.exit()
        
        # Pull in common orthologs paths
        for fname in os.listdir(os.path.join(self.project_dir, "species_data")):
            if fname.startswith("common_orthologs_{}_vs".format(self.species_name)):
                fpath = os.path.join(self.project_dir, "species_data", fname)
                if os.path.isfile(fpath):
                    self.common_orthologs_paths.append(fpath)
        

def build_species_info(project_dir):

    # Load in species information table
    species_info_file = os.path.join(project_dir, "species_data", "species_information_table.tsv")
    species_df = pd.read_csv(species_info_file, sep='\t')
    species_df = species_df[species_df["Genes with >= 1 phenotype edge"] > 0] # Only want species whith non zero phenotype information
    
    # Pull out taxon_id information (twice as two separate variables)
    ids_to_name = {sp_id:"-".join(sp_name.split(" ")) for sp_id, sp_name in zip(list(species_df["Taxon ID"]), list(species_df["Taxon label"]))}
    org_taxon_ids = copy.copy(ids_to_name)

    # Generate a dictionary with multiple keys pointing to the same species object
    species = {}
    for n in ids_to_name.keys():
        t_id = n.split(":")[1] # Without NCBITaxon: prefix
        species_obj = Species(t_id , ids_to_name[n], project_dir)
        species_obj.initiate_filepaths()
        species.update({n:species_obj})
        species.update({ids_to_name[n]:species_obj})
        species.update({n:species_obj})    # Add keys with NCBITaxon: prefix
        species.update({t_id:species_obj}) # Add keys without NCBITaxon: prefix
    
    return org_taxon_ids, species


def validate_species_arguments(species_dict, org_taxon_ids, input_taxon_id, prediction_network):
    """
    Ensures that the input taxon id is valid and can be used for downstream analysis.
    Will convert / reformat if necessary... (i.e. NCBITaxon:9606 --> 9606)
    """

    # Ensure input species id is valid
    input_taxon_id = str(input_taxon_id)
    if input_taxon_id not in species_dict:
        print('- ERROR, relevant taxon id must be supplied for taxon_id argument. Exiting...')
        sys.exit()

    # Ensure input species has necessary disease/phenotype information available
    fpath = getattr(species_dict[input_taxon_id], "gene_to_{}_path".format(prediction_network))
    if fpath == None:
        print("- ERROR, Species {} is missing gene_to_{} file... Exiting".format(input_taxon_id, prediction_network))
        sys.exit()
    
    # Return the ncbi taxon id NUMBER (i.e. 9606 instead of NCBITaxon:9606)
    return species_dict[input_taxon_id].taxon_id


###########################################################################################
### Config initialization functions for random, real, and phenolog ranking calculations ###

# Used in phenologs calculations for randomized trial versions
def initiate_random_species_comparison_configs(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    - Input species is compared against all other species within the monarch kg
    """

    # Initiates species information and ensures necessary files are present
    org_taxon_ids, species_dict = build_species_info(input_args.project_dir)
    t_id = validate_species_arguments(species_dict, org_taxon_ids, input_args.taxon_id, input_args.prediction_network)
    
    # Make all comparisons necessary across relevant species
    species_obj = species_dict[t_id]
    configs = []
    for common_orth_path in species_obj.common_orthologs_paths:

        # Our species names
        a_name = species_obj.species_name
        b_name = common_orth_path.split("_vs_")[1].split(".")[0]

        # Base level config
        apath = getattr(species_obj, "{}_to_ortholog_path".format(input_args.prediction_network))
        bpath = getattr(species_dict[b_name], "phenotype_to_ortholog_path")
        cpath = common_orth_path

        # Ensure all files exist
        if (not os.path.isfile(apath)) or (not os.path.isfile(bpath)) or (not os.path.isfile(cpath)):
            print("- ERROR, Species {} vs {} is missing common orthologs and or phenotype files...".format(a_name,
                                                                                                            b_name))
            sys.exit()
        
        # We want to output our random trial data to species specific directories.
        # Makes downstream analysis easier to deal with
        top_level_random_dir = os.path.join(input_args.project_dir, "random_trials")
        random_trial_species_dir = os.path.join(top_level_random_dir,
                                                "{}_{}_trials".format(a_name.replace("-", "_"), 
                                                                      input_args.prediction_network))
        if not os.path.isdir(random_trial_species_dir):
            os.makedirs(random_trial_species_dir, exist_ok=True)

        # Generate two configs for comparison (Species A-->B, and B-->A)
        config_a = {"species_a":a_name,
                    "species_b":b_name,
                    "species_a_phenotype_path":apath,
                    "species_b_phenotype_path":bpath,
                    "common_orthologs_path":cpath,
                    "output_directory":random_trial_species_dir}
        
        # Perform comparison in single direction
        configs.append(config_a)
    
    print("- {} relevant pairwise comparisons configurations created for {} species id...".format(len(configs), input_args.taxon_id))
    return t_id, configs
    

# Used in computing real phenologs calculations and leave xyz out validations
def initiate_phenologs_species_comparison_configs(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    - Input species is compared against all other species within the monarch kg
    """

    # Initiates species information and ensures necessary files are present
    org_taxon_ids, species_dict = build_species_info(input_args.project_dir)
    t_id = validate_species_arguments(species_dict, org_taxon_ids, input_args.taxon_id, input_args.prediction_network)
    
    # Make all comparisons necessary across relevant species
    species_obj = species_dict[t_id]
    configs = []
    for common_orth_path in species_obj.common_orthologs_paths:

        # Our species names
        a_name = species_obj.species_name
        b_name = common_orth_path.split("_vs_")[1].split(".")[0]

        # Base level config
        apath = getattr(species_obj, "{}_to_ortholog_path".format(input_args.prediction_network))
        bpath = getattr(species_dict[b_name], "phenotype_to_ortholog_path")
        cpath = common_orth_path
        ag2p_path = getattr(species_obj, "gene_to_{}_path".format(input_args.prediction_network))
        bg2p_path = getattr(species_dict[b_name], "gene_to_phenotype_path")

        # Ensure all files exist
        if (not os.path.isfile(apath)) or \
           (not os.path.isfile(bpath)) or \
           (not os.path.isfile(cpath)) or \
           (not os.path.isfile(ag2p_path)) or \
           (not os.path.isfile(bg2p_path)):
           print("- ERROR, Species {} vs {} is missing common orthologs and or phenotype files...".format(a_name, b_name))
           sys.exit()
        
        # We want to output our results data to species specific directories.
        # Makes downstream analysis easier to deal with
        top_level_res_dir = os.path.join(input_args.project_dir, "phenologs_results")
        results_species_dir = os.path.join(top_level_res_dir, 
                                            "{}_{}_results".format(a_name.replace("-", "_"), 
                                                                    input_args.prediction_network))

        if not os.path.isdir(results_species_dir):
            os.makedirs(results_species_dir, exist_ok=True)

        # Generate two configs for comparison (Species A-->B, and B-->A)
        config_a = {"species_a":a_name,
                    "species_b":b_name,
                    "species_a_phenotype_path":apath,
                    "species_b_phenotype_path":bpath,
                    
                    "species_a_g2p_path": ag2p_path,
                    "species_b_g2p_path": bg2p_path,
                    
                    "common_orthologs_path":cpath,
                    "output_directory":results_species_dir,
                    "prediction_network_type":input_args.prediction_network}
                    ##"fdr_path":check_file2} # This is the fdr table file computed from the randomized trials

        # Perform comparison in single direction
        configs.append(config_a)
            
    print("- {} relevant pairwise comparisons configurations created for {} species id...".format(len(configs),
                                                                                                  input_args.taxon_id))
    return t_id, configs


# Creates compute configs from input args for ortholog-->phenotype rankings (used in 06 step in pipeline and xvalidation)
def initiate_ortholog_to_phenotype_ranking_calculation_config(input_args):
    """
    - Attempts to ensure filepaths required for all computations are resolved before hand, so that
      calculations don't fail part way through. 
    - Creates pairwise comparison configuration data structures from
      input arguments. Either all comparisons or a select set from a comma seperated list of taxon ids
    - Input species is compared against all other species within the monarch kg
    """

    # Initiates species information and ensures necessary files are present
    org_taxon_ids, species_dict = build_species_info(input_args.project_dir)
    t_id = validate_species_arguments(species_dict, org_taxon_ids, input_args.taxon_id, input_args.prediction_network)
    species_obj = species_dict[t_id]

    # Format filenames and paths for fdr table
    sp_name = species_dict[t_id].species_name
    nformatted = sp_name.replace("-", "_")
    res_name = "{}_{}_results".format(nformatted, input_args.prediction_network)
    res_dir = os.path.join(input_args.project_dir, "phenologs_results", res_name)
    fdr_path = os.path.join(res_dir, "{}_fdr_table.tsv".format(sp_name))
    
    # Check if necessary fdr table exists and gather remaining filepaths necessary
    if not os.path.isdir(res_dir) or not os.path.isfile(fdr_path):
        print('- ERROR, relevant taxon id results not found for {}... Exiting'.format(sp_name))
        sys.exit()

    # Our species specific table filepaths from initial steps of pipeline
    p2o_path = getattr(species_obj, "{}_to_ortholog_path".format(input_args.prediction_network))
    g2o_path = getattr(species_obj, "gene_to_ortholog_path")
    g2p_path = getattr(species_obj, "gene_to_{}_path".format(input_args.prediction_network))

    # Pooled phenolog file name (significance cutoff included in name)
    sig_phenologs_outname = "{}_pooled_phenologs_fdr{}.tsv".format(sp_name, input_args.fdr)
    sig_phenologs_outpath = os.path.join(res_dir, sig_phenologs_outname)
    
    # Generate dictionary to hold our filepaths for easy processing downstream
    sp_file_info = {"project_dir":input_args.project_dir,
                    "results_dir":res_dir,
                    "taxon_id":t_id,
                    
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
    
    # If we supplied xtaxon_ids argument, we add that information here
    # Argument for including only select species during xvalidation ortholog-->phenotype distance calculations
    # (i.e. so we can observe how each species contributes to the predictions
    select_taxon_ids = {}
    try:
        for tx_id in input_args.xtaxon_ids.split(','):

            # Add taxon-name... For example 9606-->Homo-sapiens
            if tx_id in species_dict:
                select_taxon_ids.update({species_dict[tx_id].species_name:''})
            else:
                print('- ERROR, relevant taxon id must be supplied for xtaxon_ids argument. Exiting...')
                sys.exit()

    # No xtaxon_id argument supplied so we don't do anything
    except:
        dummyvar = 1
    
    sp_file_info.update({"xtaxon_ids":select_taxon_ids})
    return sp_file_info


####################################
### Helper and utility functions ###

# For multiprocessing
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


# Distance / weighting metric conversion of hypergeomatric parameters to arrays for pearson correlation calc
# If we really want to include pearson "distance" metric from the 2013 paper, then we need to compute by hand
# from the hg parameters that we generate in order to maintain the performance in terms of speed...
# This is currently very slow because it expands back out, and then uses scipy to calculate
def expand_hg_params_to_binary_pearson(hg_param_set):
    """
    Takes list hypergeomatric test parameters (c, N, m, n)
    c = overlap count
    N = total number of elements in the "bag"
    m = number of success marbles in the bag (number of 1s in array2) 
    n = number of times we draw from the bag (number of 1s in array1)
    """
    
    ## hg_param_set should be in the following format 
    ## [hg_overlap_counts, hg_world_counts, hg_b_counts, hg_a_counts]
    
    array0, array1 = [],[]
    if hg_param_set[0] > 0:
        
        # Fill in overlapping values....
        for i in range(0, hg_param_set[0]):
            array0.append(1)
            array1.append(1)
        
        c1 = hg_param_set[0]
        c2 = hg_param_set[0]
        for i in range(0, hg_param_set[1]-hg_param_set[0]):
            
            # Fill in ones for array0 # Note, last element in parameter set is the one we want to use here
            if c1 < hg_param_set[3]:
                c1 += 1
                array0.append(1)
                array1.append(0)
            
            # Fill in ones for array1
            elif c2 < hg_param_set[2]:
                c2 += 1
                array0.append(0)
                array1.append(1)
            
            else:
                array0.append(0)
                array1.append(0)
    
    return array0, array1
    #return pearsonr(array0, array1).pvalue

  
# Computes hypergeometric test from an input list of parameter sets [(c,N,m,n), ...]
def bulk_compute_hyper_geom(params):
    """
    Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
    c = number of common orthologs between phenotypes (ortholog matches)
    N = total number of orthologs shared between species
    m = number of orthologs in species B phenotype (how many desired objects are in bag)
    n = number of orthologs in species A phenotype (how many times we draw from bag)

    params argument should be in the following format 
    params = [(c,N,m,n), (c,N,m,n), ...] List of tuples where each is a set of hg parameters to compute
    """

    return {pr:float(hypergeom.pmf(*pr)) for pr in params} # Expand hg argments using the "*" character


# Expands the hyper geometric parameter sets into arrays for pearson correlation instead via scipy
def bulk_compute_pearson_from_hg_params(hg_params):
    pr_pvals, pr_coeffs = {}, {}
    for pr in hg_params:

        # Note, Could speed up dramatically by computing by hand instead of expanding and using scipy
        x, y = expand_hg_params_to_binary_pearson(pr) 
        res = pearsonr(x, y)
        pr_pvals.update({pr:res.pvalue})
        pr_coeffs.update({pr:res.correlation})
    
    return pr_pvals, pr_coeffs


# For post initial phenologs comparison calculations (last portions of the analysis pipeline) 
def load_fdr_table_to_lookup(fdr_path):
    """
    Loads fdr .tsv file to {(species a, species b, fdr):pvalue, ...} lookup datastructure
    """

    # Load fdr table for this species set of data and generate pvalue lookuptable from the computed fdr values
    fdr_df = pd.read_csv(fdr_path, sep='\t')
    fdr_lookup = {}
    for cname,pval,fdrval in zip(list(fdr_df["species_comparison"]), 
                                 list(fdr_df["pvalue_cutoff"]), 
                                 list(fdr_df["fdr"])):

        spa,spb = cname.split(',')
        key = (spa, spb, round(float(fdrval), 4))
        fdr_lookup.update({key:float(pval)})

    return fdr_lookup


# Takes _all_phenologs.tsv results files from within the results_dir, 
# pools and filters by relevant pvalue cutoff for input fdr_levle (.95 is default (i.e. 5%))
def pool_phenologs_data(results_dir, fdr_lookup_table, sig_outpath, fdr_level:.95, compress:bool=False):
    """
    Important note... The fdr lookup table is used to determin which value we need to use for each 
    _all_phenologs.tsv file that is present 
    (i.e. which fdr-pvalue cutoff we need to use for each SpeciesA-SpeciesB comparison)
    """
        
    # Assemble relevant phenologs results filepaths
    phenologs_paths = [os.path.join(results_dir, fname) for fname in os.listdir(results_dir) if fname.endswith("_all_phenologs.tsv")]
    if len(phenologs_paths) == 0:
        print("- Error, No files ending with _all_phenologs.tsv found in {}... Exiting".format(results_dir))
        sys.exit()
    
    
    # Loop through all x species phenologs files and pool together the data that passes pvalue criteria
    sig_phenologs_df = {}
    df_init_var = 0
    for rpath in phenologs_paths:
        
        
        # Generate fdr-pvalue to filter data by from fdr table file by filename scheme formatting
        # Filename should look something like this speci-es-A_vs_sp-ciesB_all_phenologs.tsv)
        spa, spb = rpath.split('/')[-1].replace("_all_phenologs.tsv", "").split("_vs_")
        fdr_key = (spa, spb, fdr_level)
        if fdr_key not in fdr_lookup_table:
            print("- ERROR, Non existant fdr lookup value... Try .95, .99, .999, .9999, etc... Exiting")
            sys.exit()
        pval_max = fdr_lookup_table[fdr_key]
        
        # Load phenologs data and filter for significant phenologs only
        phenologs_df = pd.read_csv(rpath, sep='\t')
        org_phenolog_count = len(list(phenologs_df["hg_pval"]))
        
        phenologs_df["hg_pval"] = phenologs_df["hg_pval"].astype(float)
        phenologs_df = phenologs_df[phenologs_df["hg_pval"] < pval_max]
        filtered_phenolog_count = len(list(phenologs_df["hg_pval"]))
        
        # Create / add new column to dataframe that descibes the cross species comparison that was made
        xspecies_name = ','.join([spa,spb])
        phenologs_df["X Species Comparison"] = [xspecies_name for i in range(0, filtered_phenolog_count)]
        
        # Initiate our pooled datastructure if this is the first result file we are processing
        if df_init_var == 0:
            sig_phenologs_df = phenologs_df
            df_init_var += 1
        
        # Otherwise, add our filtered phenologs to the larger pool
        else:
            sig_phenologs_df = pd.concat([sig_phenologs_df, copy.copy(phenologs_df)])
        

        # Write individual species comparison significant phenologs files
        phenologs_df = phenologs_df.sort_values(by="hg_pval")
        # Default is no gzip
        if compress == False:
            phenologs_df.to_csv(sig_outpath, sep='\t', index=False)
        else:
            phenologs_df.to_csv(sig_outpath, sep='\t', index=False, compression='gzip')


        print("- Loaded {}/{} phenolgs from {} x {} passing with pvalue <= {}".format(format(filtered_phenolog_count, ','),
                                                                                      format(org_phenolog_count, ','),
                                                                                      spa, 
                                                                                      spb, 
                                                                                      pval_max))
        
    # Write pooled results to file ###
    sig_phenologs_df = sig_phenologs_df.sort_values(by="hg_pval")

    # Default is no gzip
    if compress == False:
        sig_phenologs_df.to_csv(sig_outpath, sep='\t', index=False)
    else:
        sig_phenologs_df.to_csv(sig_outpath, sep='\t', index=False, compression='gzip')


    tot_phenologs = len(sig_phenologs_df["hg_pval"])
    print("- {} total x species phenologs pooled together...".format(format(tot_phenologs, ',')))
    print("- Data written to {}".format(sig_outpath))

    return sig_phenologs_df


# For converting / creating similarity tables between species ontologies (in phenolog space)
def convert_phenologs_to_similarity_tables(infile_path, out_dir, gzip=True):
    """
    Takes as input the filepath the species_pooled_phenologs_fdr*.tsv file and converts it to similarity tables for each species present.
    Each cross species comparison is written to out_dir as it's own separate file.
    """

    df = pd.read_csv(infile_path, sep='\t')
    comps = set(list(df["X Species Comparison"]))
    
    # Can add 1. - pvalue column here...
    df["hg_pval"] = df["hg_pval"].astype(float)

    # Pull relevant columns, and add new column for phenologs similarity
    # Make output directory if doesn't exist
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
        print("- Output directory made at {} ...".format(out_dir))
    
    # Can add 1. - pvalue column here...
    df["hg_pval"] = df["hg_pval"].astype(float)

    # Converting to similarity score via (1.0 - pvalue) has high precision requirments, and therefore an alternative method is
    # used. We take the log10 transform of our value and subtract it from 1.0. So a pvalue of 1.33e-133 --> 133.
    # We can then divide these data points by the maximum value to get everything normalized between 0.0 and 1.0 
    logsims = [1.0 - math.log(v, 10) for v in list(df["hg_pval"])]
    max_global_simscore = max(logsims)

    for c in comps:

        # Note, we format outfile path based on species name rather than ontology prefixes, 
        # because different species can use the same ontology to describe phenotypes
        spa, spb = c.split(",")
        out_name = "{}_vs_{}_sig_phenologs.tsv".format(spa, spb)
        outfile_path = os.path.join(args.out_dir, out_name)

        # Should already be sorted by pvalue... Same or other options can be applied here though 
        out_df = copy.copy(df[df["X Species Comparison"] == c])
        logsims = [1.0 - math.log(v, 10) for v in list(out_df["hg_pval"])]

        spx_max = max(logsims )
        spx_normed_global = np.asarray(logsims) / max_global_simscore
        spx_normed_spx = np.asarray(logsims) / spx_max

        # Create new df and add our two new similarity columns to our data
        default_na_col = [None for i in range(out_df.shape[0])]

        oak_df = {"subject_id":list(out_df["Species A Phenotype ID"]),
                  "subject_label":list(out_df["Species A Phenotype Name"]),
                  "subject_source":default_na_col,
                  "object_id":list(out_df["Species B Phenotype ID"]),
                  "object_label":list(out_df["Species B Phenotype Name"]),
                  "object_source":default_na_col,
                  "ancestor_id":default_na_col,
                  "ancestor_label":default_na_col,
                  "ancestor_source":default_na_col,
                  "object_information_content":default_na_col,
                  "subject_information_content":default_na_col,
                  "ancestor_information_content":default_na_col,
                  "jaccard_similarity":default_na_col,
                  "cosine_similarity":default_na_col,
                  "dice_similarity":default_na_col,
                  "phenodigm_score":default_na_col,
                  "phenolog_sim_allspecies":spx_normed_global,
                  "phenolog_sim_xspecific":spx_normed_spx}

        # Write data compressed
        if gzip == True:
            outfile_path = "{}.gz".format(outfile_path)
            pd.DataFrame(oak_df).to_csv(outfile_path, sep='\t', index=False, compression="gzip")
        
        # Write data normal
        else:
            pd.DataFrame(oak_df).to_csv(outfile_path, sep='\t', index=False)

        print("- Data written for {} --> {}".format(c, outfile_path))
    
    print("- Done!")


############################################################################
### Class's for performing phenologs calculations and validation methods ###

# For initial phenologs comparison / distance calculations
class SpeciesComparison(BaseModel):
    species_a: str # Example: Homo-sampiens, Mus-musculus, Xenopus-tropicalis...
    species_b: str
    species_a_phenotype_path: str # Example: Homo-sapiens_phenotype_to_ortholog.pkl	
    species_b_phenotype_path: str
    common_orthologs_path: str # Example: common_orthologs_Homo-sapiens_vs_Mus-musculus.tsv
    
    species_a_g2p_path: Optional[str] = None # Example Homo-sampiens_gene_to_phenotype.tsv
    species_b_g2p_path: Optional[str] = None 
    base_a_p2g: Optional[dict] = None
    base_b_p2g: Optional[dict] = None
    base_pool: Optional[list] = None
    common_ortholog_count: Optional[int] = None
    output_directory: Optional[str] = None
    
    fdr_path: Optional[str] = None
    prediction_network_type: Optional[str] = None # phenotype or disease
    leave_out_validate_set: Optional[set] = set() # Panther ortholog ids... For example set(["PTHR11482", ...])
    leave_out_validate_fdr: Optional[float] = .95 # FDR lookup value to use to lookup pvalue cutoff we need to use
    
    
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
     
    
    def subset_species_to_common_orthologs(self):
        """
        According to the 2013 implementaiton of phenologs, we want to be using the concept of orthogroup
        for orthology. We do NOT want 1:1 orthology representation, but rather collapse to protein family.
        Therefore, all gene ids (subfamilies) of the same phylogenetic tree, get collapsed to a single 
        protein family node (or otholog node). All phenotypes linking to genes beloning to the same 
        protein family will now link to the single ortholog node that is created. 

        However, this leaves the question of "frequency" of phenotype terms per protein family, meaning
        the same phenotype can map to multiple gene ids within the same protein family.
        Current implmentation collapses the frequency data to 1 for all terms. 

        - Only include genes (ortholog ids) that have >= 1 phenotype id linked to them
        - Only include genes (ortholog ids) that are common between both species
        - Hyper geometric paramters can then be computed
          - Number of common orthologs between two species is the world size parameter
          - The draw size, and number of potentially matching terms can be computed from the 
            number of orthologs associated with a particular phenotype in species a, and species b
          - The "success" parameter is simply the number of the common ortholog ids found in the underlying
            gene networks between species a phenotype, and species b phenotype
        
        - Note, lists are used to represent orthologs because this is the required input type for random sampling
        """

        # Common orths is common orthologs between the two species
        common_orths = set(self.base_pool)
        a_sub, b_sub = {}, {}

        for k,v in self.base_a_p2g.items():

            # For model cross validation. Default is include everything
            # Note, our global "commonality" parameter doesn't change. This is because is theory, the set of 
            # orthologs we are leaving out still exist between the two species, its just the edges that get removed.
            # Default is empty set so everything relevant is included
            orth_list = list(set([vv for vv in v if (vv in common_orths) and (vv not in self.leave_out_validate_set)]))
            if len(orth_list) > 0:
                a_sub.update({k:orth_list})

        for k,v in self.base_b_p2g.items():
            orth_list = list(set([vv for vv in v if vv in common_orths]))
            if len(orth_list) > 0:
                b_sub.update({k:orth_list})

        return common_orths, a_sub, b_sub


# Takes advantage of the fact that we don't actually care about the in depth nature of the
# randomized results. So we can compress and combine calculations across trials to save compute
class RandomSpeciesComparison(SpeciesComparison):
    super(SpeciesComparison)
    
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
        
        # Subset only the common orthologs for each species
        # a_sub,b_sub are in the form of {phenotype:[ortholog_id, ortholog_id], ...}
        # common_orths is set(list(self.base_pool)) # Same number of entries as self.base_pool
        common_orths, a_sub, b_sub = self.subset_species_to_common_orthologs()

        #print("- {}, {}".format(len(a_sub), len(b_sub)))
        #print("- {}".format(self.common_ortholog_count))
        #print("- {}, {}".format(max([len(v) for v in a_sub.values()]), max([len(v) for v in b_sub.values()])))

        # Generate randomized dataset (Using only common orthologs from each species that have phenotypes)
        randomized_a = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in a_sub.items()}
        randomized_b = {phen:random.sample(self.base_pool, len(orths)) for phen,orths in b_sub.items()}
        

        # Convert to compact list of lists data structure (instead of big numpy array or sparse array)
        b_matrix = [v for v in randomized_b.values()]

        # Precompute / grab partial hyper geometric paramter lookups
        b_matrix_lengths = {i:len(v) for i,v in enumerate(b_matrix)}
        hg_a_count_params = {k:len(v) for k,v in randomized_a.items()}
        hg_b_count_params = {k:len(v) for k,v in randomized_b.items()}

        # Map each unique value (ortholog id) in our data to the rows it belongs to (pre processing step)
        # Crucial piece of code in allowing fast computation of gene network overlaps between cross species phenotypes.
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
        
        # For keeping track of our hg params
        hg_a_counts = []
        hg_b_counts = []
        hg_overlap_counts = []
        hg_world_counts = []
        
        # Loop through all phenotype-ortholog associatoins for species a, and compute overlap with species b 
        for a_phen, a_orthologs in randomized_a.items():
            
            # Initialize datastructure to count the number of overlaps found for each of species b phen-orthologs
            counts = np.zeros(b_phenotype_count)
            a_ortho_count = hg_a_count_params[a_phen]
            
            # Compute commonality between each ortholog and species b orthologs.
            # We are able to compare the entiriety of species b's phenotypes using a single lookup
            # and then tallying the results. 
            for orth in a_orthologs:
                if orth in orth_to_coords:
                    counts[orth_to_coords[orth]] += 1
                t_comps += 1
            
            # This gives back an oddly shaped data array that we can flatten to get all non-zero row indices
            inds = np.argwhere(counts > 0).flatten()
            ind_count = len(inds)
            
            # World size parameter is set to the number of common orthologs between species a and b.
            # Therefore, our phenotype<-->gene networks must ONLY consist of genes that are orthologous between a and b.
            # In reality, species a phenotypeXYZ can be associated with thousands of genes (but only a handful of them
            # might be common orthologs). If we were to take all genes from a and b (not just common orthologs) then
            # we could wind up in a situation where the number of genes associated with XYZ phenotype (in species a or b)
            # is MORE than the number of orthologs those species share in common. By doing it that way, you not only
            # run into an edge case that will break that hyper geometric test, but it doesn't make sense to include
            # all genes, given our initial world size parameter is derived from only the common orothologs, and not
            # the totality of species a or species b genes. Should only be pulling from the common pool of orthologs
            # between them and the hyper geometric test should take care of the rest

            # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
            # c = number of common orthologs between phenotypes (ortholog matches)
            # N = total number of orthologs shared between species
            # m = number of orthologs in species B phenotype
            # n = number of orthologs in species A phenotype
            
            hg_overlap_counts += list(counts[inds])
            hg_world_counts += [self.common_ortholog_count] * ind_count
            hg_b_counts += [b_matrix_lengths[ind] for ind in inds]
            hg_a_counts += [a_ortho_count] * ind_count

            #processed += 1
            #if processed % 1000 == 0:
                #print("- Processed {}/10,000".format(processed))

        
        # Combine our data into a list of hg params
        tt = np.asarray([hg_overlap_counts, hg_world_counts, hg_b_counts, hg_a_counts]).T.astype(int).tolist()
        #print("- Total non zero hg params computed {}".format(format(len(tt), ',')))
        
        # Create unique set of tuples
        hg_formatted = list(map(tuple, tt))
        hg_params = set(hg_formatted)
        hg_param_counts = Counter(hg_formatted)
        #print("- Unique hg params computed {}".format(format(len(hg_params), ',')))

        return hg_params, hg_param_counts
        

    def run_randomized_comparison_trials(self, n_trials: List[int]):
        
        # Keep track of unique sets of hyper geometric tests to compute so we don't recalculate
        hg_params = set()
        data = {}
        for i, n_trial in enumerate(n_trials):

            # Generate randomized data trials
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

        # TO DO?: Is the the proper framing of this calculation / proper version of pearson we are computing... 
        pears_pvals, pears_coeffs = bulk_compute_pearson_from_hg_params(hg_params)

        
        # Map to pvalues, and write data files
        for trial_num, trial_data in data.items():
            
            # Define our outfile path and map params to pvalues 
            outfile_path = os.path.join(self.output_directory, 
                                        "{}_vs_{}_{}.tsv.gz".format(self.species_a, self.species_b, trial_num))
                                        
            # Map data back to computed pvalues, and format our trial data for easy writing to file
            pvals1 = [hg_pvals[tuple(k)] for k in trial_data]
            pvals2 = [pears_pvals[tuple(k)] for k in trial_data]
            coeffs = [pears_coeffs[tuple(k)] for k in trial_data]
            #param_data = np.asarray(trial_data).T
            
            # Reduced file size version
            param_data = np.asarray(list(trial_data.keys())).T
            occurrence = list(trial_data.values())
            
            # Write data using pandas
            pd.DataFrame({"a_ortholog_count":param_data[0],
                          "b_ortholog_count":param_data[1],
                          "overlap_count":param_data[2],
                          "common_ortholog_count":param_data[3],
                          "hg_pval":pvals1, ## Full file size version).to_csv(outfile_path, sep='\t', index=False)
                          "pearson_pval":pvals2,
                          "pearson_coeff":coeffs,
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
    

# For initial phenologs comparison / distance calculations
class PhenologsSpeciesComparison(SpeciesComparison):
    super(SpeciesComparison)
    
    
    def get_phenolog_params(self, leave_out_a=set()):
        """
        Pairwise comparisons between inter species phenotypes are made by computing the
        number of orthologs that overlap between the underlying gene networks. Note that
        ths does NOT compute pvalues, but rather gathers the set(s) of parameters used
        to calculate pvalues and full dataset of non zero overlap count parameters.
        """

        # Load initial data and reduce down to relevant sets of phenotypes / common orthologs.
        # Only phenotypes that link to one or more orthologs are included for each species
        self.load_species_ortholog_information()
        common_orths, a_sub, b_sub = self.subset_species_to_common_orthologs()
        a_phens = list(a_sub.keys())
        b_phens = list(b_sub.keys())


        # Convert to compact list data structure 
        # (instead of big numpy array or sparse array which introduce lots of uncessary overhead here)
        b_matrix = [v for v in b_sub.values()]
        b_matrix_lengths = {i:len(v) for i, v in enumerate(b_matrix)}

        hg_a_count_params = {k:len(v) for k,v in a_sub.items()}
        hg_b_count_params = {k:len(v) for k,v in b_sub.items()}

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

        # For keeping track of our hg params
        hg_a_counts = []
        hg_b_counts = []
        hg_overlap_counts = []
        hg_world_counts = []
        

        # Loop through all phenotype-ortholog associatoins for species a, and compute overlap with species b 
        species_a_phens, species_b_phens = [], []
        for a_phen, a_orthologs in a_sub.items():

            # Initialize datastructure to count the number of overlaps found for each of species b phen-orthologs
            counts = np.zeros(b_phenotype_count)
            ind_count = len(counts)
            a_ortho_count = hg_a_count_params[a_phen]

            # Compute commonality between each ortholog and species b orthologs
            for orth in a_orthologs:
                if orth in orth_to_coords:
                    counts[orth_to_coords[orth]] += 1
                t_comps += 1
            
            # This gives back an oddly shaped data array that we can flatten to get all non-zero row indices
            inds = np.argwhere(counts > 0).flatten()
            ind_count = len(inds)

            # Only tallying up non zero overlap comparisons (saves time / space)
            hg_overlap_counts += list(counts[inds].astype(int))
            hg_b_counts += [b_matrix_lengths[ind] for ind in inds]
            hg_a_counts += [a_ortho_count] * ind_count
            hg_world_counts += [self.common_ortholog_count] * ind_count
            species_a_phens += [a_phen] * ind_count
            species_b_phens += [b_phens[ind] for ind in inds]

            # This is if we want ALL data ()
            #hg_overlap_counts += list(counts.astype(int))
            #hg_b_counts += [b_matrix_lengths[ind] for ind in range(0, b_phenotype_count)]
            #hg_a_counts += [a_ortho_count] * b_phenotype_count
            #hg_world_counts += [self.common_ortholog_count] * b_phenotype_count
            #species_a_phens += [a_phen] * b_phenotype_count
            #species_b_phens += b_phens

            #processed += 1
            #if processed % 1000 == 0:
                #print("- Processed {}/10,000".format(processed)

        # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
        # c = number of common orthologs between phenotypes (ortholog matches)
        # N = total number of orthologs shared between species
        # m = number of orthologs in species B phenotype
        # n = number of orthologs in species A phenotype
        
        # Combine our data into a list of hg params
        tt = np.asarray([hg_overlap_counts, hg_world_counts, hg_b_counts, hg_a_counts]).T.astype(int).tolist()
        
        # Create unique set of tuples
        hg_formatted = list(map(tuple, tt))
        hg_params = set(hg_formatted)

        # Collate data into dictionary (easy format to work with downstream)
        return_data = {"Species A Phenotype ID":species_a_phens,
                       "Species B Phenotype ID":species_b_phens,
                       "Ortholog Count A":hg_a_counts,
                       "Ortholog Count B":hg_b_counts,
                       "Overlap Count":hg_overlap_counts,
                       "Common Ortholog Count":hg_world_counts}

        return hg_params, return_data
    
    
    def compute_cross_species_phenologs(self):
        
        # Read in fdr table to get pvalue cutoffs
        ##fdr_info = self.get_fdr_info_from_fdr_file()

        # Read in gene_to_phenotype file to get phenotype names for each species (for writing to final tables)
        species_a_g2p_df = pd.read_csv(self.species_a_g2p_path, sep='\t')
        species_b_g2p_df = pd.read_csv(self.species_b_g2p_path, sep='\t')

        # Create mapping tables of phenotype:id --> phenotype name
        # Note, we need to specify what type of predictions we are trying to make (phenotype or disease)
        akey = self.prediction_network_type
        ap2name = {k:v for k,v in zip(list(species_a_g2p_df[akey]), list(species_a_g2p_df["{}_name".format(akey)]))}
        bp2name = {k:v for k,v in zip(list(species_b_g2p_df["phenotype"]), list(species_b_g2p_df["phenotype_name"]))}

        # Compute parameters for hyper geometric tests
        hg_params, comparison_data = self.get_phenolog_params()
        print("- {} vs. {} -- hg parameters computed {}...".format(self.species_a,
                                                                   self.species_b,
                                                                   format(len(hg_params), ',')))
        
        # Compute hyper geometric tests
        hg_pvals = bulk_compute_hyper_geom(hg_params)

        # TO DO?: Is the the proper framing of this calculation / proper version of pearson we are computing... 
        #pears_pvals, pears_coeffs = bulk_compute_pearson_from_hg_params(hg_params)
        
        pval_hg_col, sig_col = [], []
        pval_pe_col, coeff_pe_col = [], []
        apname_col, bpname_col = [], []

        # Map each row's set of hg params back to the repsective pvalue (Note -- c,N,m,n param ordering matters)
        for a_phid, b_phid, c_param, N_param, m_param, n_param in zip(comparison_data["Species A Phenotype ID"],
                                                                      comparison_data["Species B Phenotype ID"],
                                                                      comparison_data["Overlap Count"], 
                                                                      comparison_data["Common Ortholog Count"],
                                                                      comparison_data["Ortholog Count B"],
                                                                      comparison_data["Ortholog Count A"]):
            
            # Note that ordering here matters
            key = tuple((c_param, N_param, m_param, n_param))
            pval_hg = hg_pvals.get(key, 1.0)
            #pval_pe = pears_pvals[key]
            #coeff_pe = pears_coeffs[key]

            # Map phenotype id to human readable name
            phena_name = ap2name[a_phid]
            phenb_name = bp2name[b_phid]

            # Add info to new output data "columns"
            pval_hg_col.append(pval_hg)
            #pval_pe_col.append(pval_pe)
            #coeff_pe_col.append(coeff_pe)

            apname_col.append(phena_name)
            bpname_col.append(phenb_name)

        # Add pvalues and convert to data frame
        comparison_data.update({"hg_pval":pval_hg_col})
        #comparison_data.update({"pearson_pval":pval_pe_col})
        #comparison_data.update({"pearson_coeff":coeff_pe_col})

        comparison_data.update({"Species A Phenotype Name":apname_col})
        comparison_data.update({"Species B Phenotype Name":bpname_col})
        comparison_data = pd.DataFrame(comparison_data)

        fname = "{}_vs_{}_all_phenologs.tsv".format(self.species_a, self.species_b)
        outpath_name = os.path.join(self.output_directory, fname)
        comparison_data.sort_values("hg_pval").to_csv(outpath_name, sep='\t', index=False)
        
        # Write select pvalue / fdr cuttoff phenologs
        #for pval_cutoff,f_ext in fdr_info.items():

            # Format output filenames
        #    fname = "{}_vs_{}_phenologs_{}.tsv".format(self.species_a, self.species_b, f_ext)
        #    outpath_name = os.path.join(self.output_directory, fname)

            # Subset df and write data

        #    comparison_data[comparison_data["P-Value"] <= pval_cutoff].sort_values("P-Value").to_csv(outpath_name, sep='\t', index=False)
            #comparison_data[comparison_data["P-Value"] <= pval_cutoff].to_csv(outpath_name, sep='\t', index=False)
        #    print("- {} written...".format(outpath_name))
        
        # TO DO: Do we need this? It's a bit excessive in terms of filesize...
        # Write all data
        ##outpath_name = os.path.join(out_directory, "{}_vs_{}_fulldata.tsv".format(self.species_a, self.species_b))                                                                   
        ##comparison_data.to_csv(outpath_name)
        return 


# For computing gene-->phenotype ranking "matrices" (last portions of the analysis pipeline)
# We make this an object, so we can multiprocess easier leave xyz out validation results
class OrthologToPhenotypeCalculations(BaseModel):
    
    # Arguments passed in from command line
    project_dir:str
    results_dir:str
    prediction_network:str # disease or phenotype
    taxon_id:str
    
    fdr_path:str
    fdr:float=.95
    kneighbs:int=10
    rank_metric:str
        
    sig_phenologs_path:str
    phen_to_orth_path:str
    gene_to_orth_path:str
    gene_to_phen_path:str
    species_name:str

    # For filtering out undesired species
    xtaxon_ids: Optional[Dict] = None
    make_new_results_dir: Optional[bool] = False
    
    
    def build_ortholog_to_phenotype_data(self):
    
        # Load phenotype (or disease) to ortholog datastructure {phenotype_id:[ortholog_id, ...]}
        p2o = pickle.load(open(self.phen_to_orth_path, 'rb'))

        # Identify unique set of orthologs available for this species from loaded data
        unique_orths = {k:'' for orths in p2o.values() for k in orths}
        print("- {} unique phenotypes found from {}".format(format(len(unique_orths), ','), 
                                                            self.phen_to_orth_path))

        # Build "matrix" linking each ortholog (protein family in this case) to the avialable set of phenotypes
        #o2p_dists = {g:{p:1. for p in p2o.keys()} for g in unique_orths}
        o2p_dists = {g:{} for g in unique_orths} 

        return o2p_dists
    
    
    def build_species_phenotype_to_orths_data(self, prediction_network):
    
        species_data_dir = os.path.join(self.project_dir, "species_data")
        species_data = {}
        for fname in os.listdir(species_data_dir):

            suffi = "_{}_to_ortholog.pkl".format(prediction_network)
            if not fname.endswith(suffi):
                continue

            # Load phenotype/disease to ortholog data, and convert to set() instead of list for easy downstream processing
            # Note, we are collapsing repeat orthologs for any given phenotype 
            # This is one of the benefits of the orthogroup abstraction, although using "frequency" information 
            # (i.e how many times a protein family/ortho_id shows up for this phenotype) may still be possible...
            spname = fname.replace(suffi, "")
            p2o = pickle.load(open(os.path.join(species_data_dir, fname), 'rb'))
            for k in p2o:
                p2o[k] = set(p2o[k])

            species_data.update({spname:p2o})
            print("- {} read into memory...".format(fname))

        # Give a warning if nothing is found...
        if len(species_data) == 0:
            print("- Error, no relevant {}_to_ortholog.pkl files found in {}".format(prediction_network, species_data_dir))

        return species_data
    

    def compute_ortholog_phenotype_distances(self, sig_phens_path: Optional[str] = None):
        
        allowed_dists = {"hg":'hyper_geometric', "nb":"naive_bayes"}
        if self.rank_metric not in allowed_dists:
            print("- ERROR, hg or nb are allowed for rank_metric argument not {}".format(self.rank_metric))
            sys.exit()


        # Load all species phenotype-->ortholog files
        p2o_species = self.build_species_phenotype_to_orths_data("phenotype")

        # Compute gene-->phenotype rankings ###
        # Precomputes "gene" to "phenotype" dictionary {gene:{phenotype:distance}}, that we can fill in later
        # Equivilant to adjacency matrix, but in dictionary form where we key each row, and
        # only generate "column" values as the come up
        o2p_dists = self.build_ortholog_to_phenotype_data()

        # Read in significant phenologs table (presumably pooled from cross species data, but can be anything)
        # Can pass in phenolog file here to get all filepaths from base level config,
        # or default is to use the filepath that is passed when initializing the object
        if not sig_phens_path:
            sig_phens_path = self.sig_phenologs_path

        if sig_phens_path.endswith(".gz"):
            sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t', compression="gzip")
        else:
            sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t')

        # Build "k-nearest neighbor" data structure for each "phenotype" we want to assign gene rankings to
        # First, map row indices to each phenotype
        ind = 0
        p2_phens = {}

        xspecies_filtered = 0
        #Species A Phenotype ID  Species B Phenotype ID  Ortholog Count A        Ortholog Count B        Overlap Count   Common Ortholog Count   hg_pval Species A Phenotype Name        Species B Phenotype Name             X Species Comparison
        for phen_id, phen_dist_val, xcomp_name in zip(list(sig_phenologs_df["Species A Phenotype ID"]), 
                                                      list(sig_phenologs_df["hg_pval"]),
                                                      list(sig_phenologs_df["X Species Comparison"])):

            # Can filter out xyz species id here for comparisons of how each species affects performance
            spaid, spbid = xcomp_name.split(",")
            if self.xtaxon_ids:
                if spbid not in self.xtaxon_ids:
                    ind += 1 # Still need to update index even though we don't pull this data in
                    xspecies_filtered += 1
                    continue

            if phen_id not in p2_phens:
                p2_phens.update({phen_id:[]})
            p2_phens[phen_id].append(ind)
            ind += 1
        
        # Now make mini dataframes for each "phenotype" identifier and sort by best "distance" 
        # for whichever metric is chosen
        zero_sig_phens = {}
        for k in p2_phens:
            p2_phens[k] = copy.copy(sig_phenologs_df.iloc[p2_phens[k]])
            p2_phens[k] = p2_phens[k].sort_values(by="hg_pval")

            # Ensure data is of proper type for calculations
            p2_phens[k]["Overlap Count"] = p2_phens[k]["Overlap Count"].astype(float)
            p2_phens[k]["Ortholog Count B"] = p2_phens[k]["Ortholog Count B"].astype(float)

            # Now compute from our subsetted data, the cumulitive probability from 
            # "naive bayes" method, hg method?, or other?
            knear = 1
            prob_val_naiv_bayes = 1.
            hg_c, hg_N, hg_m, hg_n = 0, 0, 0, 0

            guilty_orths = set()
            for dist_val, overlap_val, comm_val, ortho_a_val, ortho_b_val, phen_id, sp_comp in zip(p2_phens[k]["hg_pval"], 
                                                                                                   p2_phens[k]["Overlap Count"], 
                                                                                                   p2_phens[k]["Common Ortholog Count"],
                                                                                                   p2_phens[k]["Ortholog Count A"],
                                                                                                   p2_phens[k]["Ortholog Count B"],
                                                                                                   p2_phens[k]["Species B Phenotype ID"],
                                                                                                   p2_phens[k]["X Species Comparison"]):

                # Compute probability via "naive bayes" method
                prob_val_naiv_bayes *= ( 1. - ((overlap_val/ortho_b_val)*(1. - dist_val)) )

                # Compute probability via second round of hyper geometric...
                hg_c += overlap_val
                hg_N += comm_val
                hg_m += ortho_b_val
                hg_n += ortho_a_val


                # Keep track of "implicated" genes/orthologs by leveraging lookup tables made earlier
                spa_name, spb_name = sp_comp.split(",")
                guilty_orths = guilty_orths | p2o_species[spb_name][phen_id]
                knear += 1

                if knear > self.kneighbs:
                    break

            # Note, for naive bayes method... we are computing the probability of at least one association NOT occuring 
            # within the top knearest neighbors. Then all gene orthologs associated get weighted equally  

            # Edge case for no significant orthologs associated with phenotype
            if len(guilty_orths) == 0:
                zero_sig_phens.update({k:''})
                continue


            
            # Precompute / determine what value we need to update
            # TO DO: Currently, if 100 or more neighbors are used in the calculation, the pvalues get rounded to zero
            # Because we are ranking data, we can divide by a factor of 10 so our rankings don't collapse 
            # Is this the best way to do this other than altering the methodology? 
            # It seems to perform better than the naive bayes for larger k but worse for very small k...
            # but improvements could be made
            if self.rank_metric == "nb": # Naive bayes scheme
                metric_val = prob_val_naiv_bayes
            elif self.rank_metric == "hg": # Hypergeometric scheme
                # TO DO: Better cuttoff (precompute the largest "world size" paramter
                # hypergeomtric test that doesn't collapse to 1.0 or 0.0 )
                if self.kneighbs >= 100: 
                    metric_val = float(hypergeom.pmf(int(hg_c/10.), 
                                                     int(hg_N/10.), 
                                                     int(hg_m/10.), 
                                                     int(hg_n/10.)))
                else:
                    metric_val = float(hypergeom.pmf(hg_c, hg_N, hg_m, hg_n))
            
            # Now fill in gene-phenotype "matrix"
            for gorth in guilty_orths:
                if gorth not in o2p_dists: # Not a common ortholog between the two species so we can't make a statement
                    continue
                
                # Fill in data as we need to
                if k not in o2p_dists[gorth]:
                    o2p_dists[gorth].update({k:metric_val})
        
        # Write out data here (current is .pkl dictionary, tor read in later, but might be nice to have
        # more human readable form... leave xyz out strategy also produces a lot of data so need small file sizes)
        sgpp = Path(sig_phens_path)
        n1, n2 = sgpp.name.split("_pooled_phenologs_") # Splits our filename into two parts we can combine
        
        # We need to deal with our xtaxons ids here so we can delineate output data
        if self.xtaxon_ids:
            xtids_formatted = "_".join(list(self.xtaxon_ids.keys()))
            outdir_name = "ortholog_to_{}_{}_{}kneighbs_{}_{}".format(self.prediction_network, 
                                                                      xtids_formatted, 
                                                                      self.kneighbs, 
                                                                      self.rank_metric,
                                                                      self.fdr)

        else:
            outdir_name = "ortholog_to_{}_{}kneighbs_{}_{}".format(self.prediction_network, 
                                                                   self.kneighbs, 
                                                                   self.rank_metric,
                                                                   self.fdr)


        # Write data to parent results directory (default mode)
        if self.make_new_results_dir == False:
            outfile_path = os.path.join(sgpp.parent, "{}_{}.pkl".format(n1, outdir_name))
        
        # Otherwise, make new directory within parent resutls directory
        elif self.make_new_results_dir == True:
            outdir_path = os.path.join(sgpp.parent, outdir_name)
            outfile_path = os.path.join(outdir_path, "{}_{}.pkl".format(n1, outdir_name))
            if not os.path.isdir(outdir_path):
                os.makedirs(outdir_path, exist_ok=True)
        
        # Write data to designated path/directory
        pickle.dump(o2p_dists, open(outfile_path, 'wb'))

        #return o2p_dists
        #print(format(sum([len(v) for k,v in o2p_dists.items()]), ','), sig_phens_path)
        ##collapsed_d{orth_id:{phen_id:dist_val for phen_id,dist_val  in op_dists[orth_id].items() if dist_val < 1.} for orth_id in op_dists}
        return
    

    # Experimental...
    def compute_ortholog_phenotype_distances_hg2(self, sig_phens_path: Optional[str] = None):
        
        allowed_dists = {"hg":'hyper_geometric', "nb":"naive_bayes"}
        if self.rank_metric not in allowed_dists:
            print("- ERROR, hg or nb are allowed for rank_metric argument not {}".format(self.rank_metric))
            sys.exit()


        # Load all species phenotype-->ortholog files
        p2o_species = self.build_species_phenotype_to_orths_data("phenotype")

        # Compute gene-->phenotype rankings ###
        # Precomputes "gene" to "phenotype" dictionary {gene:{phenotype:distance}}, that we can fill in later
        # Equivilant to adjacency matrix, but in dictionary form where we key each row, and
        # only generate "column" values as the come up
        o2p_dists = self.build_ortholog_to_phenotype_data()

        # Read in significant phenologs table (presumably pooled from cross species data, but can be anything)
        # Can pass in phenolog file here to get all filepaths from base level config,
        # or default is to use the filepath that is passed when initializing the object
        if not sig_phens_path:
            sig_phens_path = self.sig_phenologs_path

        if sig_phens_path.endswith(".gz"):
            sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t', compression="gzip")
        else:
            sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t')

        # Build "k-nearest neighbor" data structure for each "phenotype" we want to assign gene rankings to
        # First, map row indices to each phenotype
        ind = 0
        p2_phens = {}

        xspecies_filtered = 0
        #Species A Phenotype ID  Species B Phenotype ID  Ortholog Count A        Ortholog Count B        Overlap Count   Common Ortholog Count   hg_pval Species A Phenotype Name        Species B Phenotype Name             X Species Comparison
        for phen_id, phen_dist_val, xcomp_name in zip(list(sig_phenologs_df["Species A Phenotype ID"]), 
                                                      list(sig_phenologs_df["hg_pval"]),
                                                      list(sig_phenologs_df["X Species Comparison"])):

            # Can filter out xyz species id here for comparisons of how each species affects performance
            spaid, spbid = xcomp_name.split(",")
            if self.xtaxon_ids:
                if spbid not in self.xtaxon_ids:
                    ind += 1 # Still need to update index even though we don't pull this data in
                    xspecies_filtered += 1
                    continue

            if phen_id not in p2_phens:
                p2_phens.update({phen_id:[]})
            p2_phens[phen_id].append(ind)
            ind += 1
        
        # Now make mini dataframes for each "phenotype" identifier and sort by best "distance" 
        # for whichever metric is chosen
        zero_sig_phens = {}
        for k in p2_phens:
            p2_phens[k] = copy.copy(sig_phenologs_df.iloc[p2_phens[k]])
            p2_phens[k] = p2_phens[k].sort_values(by="hg_pval")

            # Ensure data is of proper type for calculations
            p2_phens[k]["Overlap Count"] = p2_phens[k]["Overlap Count"].astype(float)
            p2_phens[k]["Ortholog Count B"] = p2_phens[k]["Ortholog Count B"].astype(float)

            # Now compute from our subsetted data, the cumulitive probability from 
            # "naive bayes" method, hg method?, or other?
            knear = 1
            prob_val_naiv_bayes = 1.
            hg_c, hg_N, hg_m, hg_n = 0, 0, 0, 0

            guilty_orths = Counter()
            for dist_val, overlap_val, comm_val, ortho_a_val, ortho_b_val, phen_id, sp_comp in zip(p2_phens[k]["hg_pval"], 
                                                                                                   p2_phens[k]["Overlap Count"], 
                                                                                                   p2_phens[k]["Common Ortholog Count"],
                                                                                                   p2_phens[k]["Ortholog Count A"],
                                                                                                   p2_phens[k]["Ortholog Count B"],
                                                                                                   p2_phens[k]["Species B Phenotype ID"],
                                                                                                   p2_phens[k]["X Species Comparison"]):

                # Compute probability via "naive bayes" method
                prob_val_naiv_bayes *= ( 1. - ((overlap_val/ortho_b_val)*(1. - dist_val)) )

                # Compute probability via second round of hyper geometric...
                hg_c += overlap_val
                hg_N += comm_val
                hg_m += ortho_b_val
                hg_n += ortho_a_val


                # Keep track of "implicated" genes/orthologs by leveraging lookup tables made earlier
                spa_name, spb_name = sp_comp.split(",")
                for guilty_orth in p2o_species[spb_name][phen_id]:
                    guilty_orths[guilty_orth] += 1
                knear += 1

                if knear > self.kneighbs:
                    break

            # Note, for naive bayes method... we are computing the probability of at least one association NOT occuring 
            # within the top knearest neighbors. Then all gene orthologs associated get weighted equally  

            # Edge case for no significant orthologs associated with phenotype
            if len(guilty_orths) == 0:
                zero_sig_phens.update({k:''})
                continue


            
            # Precompute / determine what value we need to update
            # TO DO: Currently, if 100 or more neighbors are used in the calculation, the pvalues get rounded to zero
            # Because we are ranking data, we can divide by a factor of 10 so our rankings don't collapse 
            # Is this the best way to do this other than altering the methodology? 
            # It seems to perform better than the naive bayes for larger k but worse for very small k...
            # but improvements could be made
            if self.rank_metric == "nb": # Naive bayes scheme
                metric_val = prob_val_naiv_bayes
            elif self.rank_metric == "hg": # Hypergeometric scheme
                # TO DO: Better cuttoff (precompute the largest "world size" paramter
                # hypergeomtric test that doesn't collapse to 1.0 or 0.0 )
                if self.kneighbs >= 100: 
                    metric_val = float(hypergeom.pmf(int(hg_c/10.), 
                                                     int(hg_N/10.), 
                                                     int(hg_m/10.), 
                                                     int(hg_n/10.)))
                else:
                    metric_val = float(hypergeom.pmf(hg_c, hg_N, hg_m, hg_n))
            
            # Now fill in gene-phenotype "matrix"
            for gorth in guilty_orths:
                if gorth not in o2p_dists: # Not a common ortholog between the two species so we can't make a statement
                    continue
                
                # Fill in data as we need to
                if k not in o2p_dists[gorth]:
                    o2p_dists[gorth].update({k:metric_val})
        
        # Write out data here (current is .pkl dictionary, tor read in later, but might be nice to have
        # more human readable form... leave xyz out strategy also produces a lot of data so need small file sizes)
        sgpp = Path(sig_phens_path)
        n1, n2 = sgpp.name.split("_pooled_phenologs_") # Splits our filename into two parts we can combine
        
        # We need to deal with our xtaxons ids here so we can delineate output data
        if self.xtaxon_ids:
            xtids_formatted = "_".join(list(self.xtaxon_ids.keys()))
            outdir_name = "ortholog_to_{}_{}_{}kneighbs_{}_{}".format(self.prediction_network, 
                                                                      xtids_formatted, 
                                                                      self.kneighbs, 
                                                                      self.rank_metric,
                                                                      self.fdr)

        else:
            outdir_name = "ortholog_to_{}_{}kneighbs_{}_{}".format(self.prediction_network, 
                                                                   self.kneighbs, 
                                                                   self.rank_metric,
                                                                   self.fdr)


        # Write data to parent results directory (default mode)
        if self.make_new_results_dir == False:
            outfile_path = os.path.join(sgpp.parent, "{}_{}.pkl".format(n1, outdir_name))
        
        # Otherwise, make new directory within parent resutls directory
        elif self.make_new_results_dir == True:
            outdir_path = os.path.join(sgpp.parent, outdir_name)
            outfile_path = os.path.join(outdir_path, "{}_{}.pkl".format(n1, outdir_name))
            if not os.path.isdir(outdir_path):
                os.makedirs(outdir_path, exist_ok=True)
        
        # Write data to designated path/directory
        pickle.dump(o2p_dists, open(outfile_path, 'wb'))

        #return o2p_dists
        #print(format(sum([len(v) for k,v in o2p_dists.items()]), ','), sig_phens_path)
        ##collapsed_d{orth_id:{phen_id:dist_val for phen_id,dist_val  in op_dists[orth_id].items() if dist_val < 1.} for orth_id in op_dists}
        return
    

    ### Experimental...
    def batch_compute_ortholog_phenotype_distances(self, sig_phens_paths: list):
        
        #####################################################
        ### Load as much information up front as possible ###

        # FDR tsv table to dictionary lookup {(species a, species b, fdr):pvalue, ...} 
        fdr_lookup = load_fdr_table_to_lookup(self.fdr_path)

        # Loads "prediction" species relevent network that we want predictions for (phenotype or disease)
        sp_p2o = self.build_species_phenotype_to_orths_data(self.prediction_network)
                                                       
        # Load all species phenotype-->ortholog files
        p2o_species = self.build_species_phenotype_to_orths_data("phenotype")

        ##################################################
        ### Compute results for each filepath provided ###

        for sig_phens_path in sig_phens_paths:
            if sig_phens_path.endswith(".gz"):
                sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t', compression="gzip")
            else:
                sig_phenologs_df = pd.read_csv(sig_phens_path, sep='\t')
            
            # Compute gene-->phenotype rankings ###
            # Precomputes "gene" to "phenotype" dictionary {gene:{phenotype:distance}}, that we can fill in later
            # Equivilant to adjacency matrix, but in dictionary form...
            o2p_dists = self.build_ortholog_to_phenotype_data()

            # Build "k-nearest neighbor" data structure for each "phenotype" we want to assign gene rankings to
            # First, map row indices to each phenotype
            ind = 0
            p2_phens = {}

            #Species A Phenotype ID  Species B Phenotype ID  Ortholog Count A        Ortholog Count B        Overlap Count   Common Ortholog Count   hg_pval Species A Phenotype Name        Species B Phenotype Name             X Species Comparison
            for phen_id, phen_dist_val in zip(list(sig_phenologs_df["Species A Phenotype ID"]), 
                                              list(sig_phenologs_df["hg_pval"])):

                if phen_id not in p2_phens:
                    p2_phens.update({phen_id:[]})
                p2_phens[phen_id].append(ind)
                ind += 1
            
            
            # Now make mini dataframes for each "phenotype" identifier and sort by best "distance" 
            # for whichever metric is chosen
            zero_sig_phens = {}
            for k in p2_phens:
                p2_phens[k] = copy.copy(sig_phenologs_df.iloc[p2_phens[k]])
                p2_phens[k] = p2_phens[k].sort_values(by="hg_pval")

                # Ensure data is of proper type for calculations
                p2_phens[k]["Overlap Count"] = p2_phens[k]["Overlap Count"].astype(float)
                p2_phens[k]["Ortholog Count B"] = p2_phens[k]["Ortholog Count B"].astype(float)

                # Now compute from our subsetted data, the cumulitive probability from 
                # "naive bayes" method, hg method?, or other?
                knear = 1
                prob_val_naiv_bayes = 1.
                hg_c, hg_N, hg_m, hg_n = 0, 0, 0, 0

                guilty_orths = set()
                for dist_val, overlap_val, comm_val, ortho_a_val, ortho_b_val, phen_id, sp_comp in zip(p2_phens[k]["hg_pval"], 
                                                                                                    p2_phens[k]["Overlap Count"], 
                                                                                                    p2_phens[k]["Common Ortholog Count"],
                                                                                                    p2_phens[k]["Ortholog Count A"],
                                                                                                    p2_phens[k]["Ortholog Count B"],
                                                                                                    p2_phens[k]["Species B Phenotype ID"],
                                                                                                    p2_phens[k]["X Species Comparison"]):

                    # Compute probability via "naive bayes" method
                    prob_val_naiv_bayes *= ( 1. - ((overlap_val/ortho_b_val)*(1. - dist_val)) )

                    # Compute probability via second round of hyper geometric...
                    hg_c += overlap_val
                    hg_N += comm_val
                    hg_m += ortho_b_val
                    hg_n += ortho_a_val


                    # Keep track of "implicated" genes/orthologs by leveraging lookup tables made earlier
                    spa_name, spb_name = sp_comp.split(",")
                    guilty_orths = guilty_orths | p2o_species[spb_name][phen_id]
                    knear += 1

                    if knear > self.kneighbs:
                        break

                # Note, for naive bayes method... we are computing the probability of at least one association NOT occuring 
                # within the top knearest neighbors. Then all gene orthologs associated get weighted equally  

                # Edge case for no significant orthologs associated with phenotype
                if len(guilty_orths) == 0:
                    zero_sig_phens.update({k:''})
                    continue

                # Now fill in gene-phenotype "matrix"
                for gorth in guilty_orths:
                    if gorth not in o2p_dists: # Not a common ortholog between the two species so we can't make a statement
                        continue

                    o2p_dists[gorth][k] = prob_val_naiv_bayes
                    #o2p_dists[gorth][k] = float(hypergeom.pmf(hg_c, hg_N, hg_m, hg_n))
            
            #return o2p_dists
        return