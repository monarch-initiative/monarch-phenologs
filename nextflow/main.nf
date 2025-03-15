#!/usr/bin/env nextflow

include { get_phenologs_env } from './modules/01_get_env.nf'
include { get_phenologs_data } from './modules/02_get_data.nf'
include { compute_fdr_data } from './modules/03_compute_random_trials.nf'
include { compute_real_phenolog_data } from './modules/04_compute_phenologs.nf'
include { compute_fdr_info } from './modules/05_compute_fdr.nf'
include { compute_ortholog_rank_calcs } from './modules/06_compute_ortholog_rankings.nf'
include { leave_one_out_calculations } from './modules/leave_one_out_calculations.nf'
include { leave_one_out_ortholog_rank_calcs } from './modules/leave_one_out_rankings.nf'


// Phenologs calculation parameters 
params.n_random_trials = 2
params.cpu_cores = 10
params.taxon_id = 9606
params.prd = "disease"

// Phenologs ortholog to phenotype ranking params
params.fdr = .95
params.kneighbs = 10
params.rank_metric = 'nb'

// Xvalidation options
params.xvalidate_calc = false
params.xvalidate_rank = false


workflow {

    // Get the monarch-phenologs environment
    get_phenologs_env()
    get_phenologs_data(get_phenologs_env.out)

    // Compute random trials
    compute_fdr_data(get_phenologs_env.out, 
                     get_phenologs_data.out,
                     params)
    
    // Compute real phenolog data
    compute_real_phenolog_data(get_phenologs_env.out, 
                               get_phenologs_data.out,
                               params)
    
    // Compute FDR info for significance cut
    compute_fdr_info(get_phenologs_env.out, 
                     get_phenologs_data.out,
                     compute_fdr_data.out.comp_sig_trials,  // Completion signal for random trials
                     compute_real_phenolog_data.out.comp_sig_real,
                     params)  // Completion signal for real phenolog calcs

    // Compute ortholog to phenotype rankings
    compute_ortholog_rank_calcs(get_phenologs_env.out, 
                                get_phenologs_data.out,
                                compute_fdr_info.out.comp_sig_fdr,
                                params)  // Completion signal for FDR info 
    
    // Leave one out cross validation
    if (params.xvalidate_calc) {
        leave_one_out_calculations(get_phenologs_env.out, 
                                   get_phenologs_data.out,
                                   compute_fdr_info.out.comp_sig_fdr,  // Completion signal for FDR info
                                   params)
                                   }  
    
    // Leave one out cross validation assessment
    if (params.xvalidate_rank) {
        leave_one_out_ortholog_rank_calcs(get_phenologs_env.out, 
                                        get_phenologs_data.out,
                                        leave_one_out_calculations.out.comp_sig_xval, // Completion signal for xval calcs
                                        params)
    }

}