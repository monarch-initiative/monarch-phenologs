# General imports
import os
import sys
import argparse
import copy
import math
import pandas as pd
import numpy as np

# Custom imports
from phenologs_utils import (build_species_info, validate_species_arguments)


if __name__ == '__main__':
    ################
	## ARG PARSE ###
    def parse_input_command():
        parser = argparse.ArgumentParser(description='Computes final phenolog calculations between species. Tables \
                                                      are written for each species comparison based on fdr levels.')

        parser.add_argument("-p","--project_dir", help="Top most project directory", required=False, type=str, default=None)
        parser.add_argument("-taxon_id", help='Specicies specific taxon id or "all" are allowed', required=True, type=str)
        parser.add_argument("-prd", "--prediction_network", help="phenotype or disease (which type of network to use for base species comparisons)", required=True, default="phenotype")
        parser.add_argument("-fdr", help="One minus the false discovery rate.. .95 is default", required=True, type=float, default=.95)
        return parser.parse_args()

    args = parse_input_command()
    ############################

    # Initiates species information and ensures necessary files are present
    org_taxon_ids, species_dict = build_species_info(args.project_dir)
    t_id = validate_species_arguments(species_dict, org_taxon_ids, args.taxon_id, args.prediction_network)
    species_obj = species_dict[t_id]

    ### File name should be formated something like this... Homo-sapiens_pooled_phenologs_fdr0.95.tsv
    fname = "{}_pooled_phenologs_fdr{}.tsv".format(species_obj.species_name, args.fdr)
    fdir = os.path.join(args.project_dir, "phenologs_results", "{}_{}_results".format(species_obj.species_name.replace("-", "_"), 
                                                                                       args.prediction_network))

    # Input file and output directory
    fpath = os.path.join(fdir, fname)
    outdir = os.path.join(fdir, "similarity_tables_fdr{}".format(args.fdr))
    if not os.path.isfile(fpath):
        print("- ERROR, {} not found. Exiting...".format(fpath))
        sys.exit()
    os.makedirs(outdir, exist_ok=True)

    # Read relevant pooled datafile into memory
    df = pd.read_csv(fpath, sep='\t')
    comps = set(list(df["X Species Comparison"]))
    
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
        out_name_xall = "{}_vs_{}_allspecies.semsim.tsv.gz".format(spa, spb)
        out_name_xspecific = "{}_vs_{}_xspecific.semsim.tsv.gz".format(spa, spb)

        outfile_path_xall = os.path.join(outdir, out_name_xall)
        outfile_path_xspecific = os.path.join(outdir, out_name_xspecific)

        # Should already be sorted by pvalue... Same or other options can be applied here though 
        out_df = copy.copy(df[df["X Species Comparison"] == c])
        logsims = [1.0 - math.log(v, 10) for v in list(out_df["hg_pval"])]

        spx_max = max(logsims )
        spx_normed_global = np.asarray(logsims) / max_global_simscore
        spx_normed_spx = np.asarray(logsims) / spx_max

        # Create new df and add our two new similarity columns to our data
        default_na_col = [None for i in range(out_df.shape[0])]

        df_xall = {"subject_id":list(out_df["Species A Phenotype ID"]),
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
                  "phenodigm_score":spx_normed_global} # Note, the name of this column is for downstream application of exomiser

        df_xspecific = copy.copy(df_xall)
        df_xspecific["phenodigm_score"] = spx_normed_spx

        # Write data compressed
        pd.DataFrame(df_xall).to_csv(outfile_path_xall, sep='\t', index=False, compression="gzip")
        pd.DataFrame(df_xspecific).to_csv(outfile_path_xspecific, sep='\t', index=False, compression="gzip")
        print("- Data written for {}".format(c))
    
    print("- Done!")