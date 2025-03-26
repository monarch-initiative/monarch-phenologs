#!/usr/bin/env nextflow

include { get_phenologs_env } from './modules/01_get_env.nf'
include { get_phenologs_data } from './modules/02_get_data.nf'
include { compute_fdr_data } from './modules/03_compute_random_trials.nf'
include { compute_real_phenolog_data } from './modules/04_compute_phenologs.nf'
include { compute_fdr_info } from './modules/05_compute_fdr.nf'
include { compute_ortholog_rank_calcs } from './modules/06_compute_ortholog_rankings.nf'
include { convert_to_sim_tables } from './modules/to_similarity_tables.nf'
include { leave_one_out_calculations } from './modules/leave_one_out_calculations.nf'
include { leave_one_out_ortholog_rank_calcs } from './modules/leave_one_out_rankings.nf'


// Phenologs calculation parameters 
params.n_random_trials = 20
params.cpu_cores = 10 // Not actuall used (any more... config takes care of this)
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

    Channel.value(params.taxon_id).set{ taxon_id }
    Channel.value(params.prd).set{ prd }
    Channel.value(params.n_random_trials).set{ n_random_trials }
    Channel.value(params.fdr).set{ fdr }
    Channel.value(params.kneighbs).set{ kneighbs }
    Channel.value(params.rank_metric).set{ rank_metric }

    // Get the monarch-phenologs environment
    get_phenologs_env()
    get_phenologs_data(get_phenologs_env.out.env_path)

    // Compute random trials
    compute_fdr_data(get_phenologs_env.out.env_path, 
                     get_phenologs_data.out.project_path,
                     taxon_id,
                     prd,
                     n_random_trials)
    
    // Compute real phenolog data
    compute_real_phenolog_data(get_phenologs_env.out.env_path, 
                               get_phenologs_data.out.project_path,
                               taxon_id,
                               prd)

    // Compute FDR info for significance cut
    compute_fdr_info(get_phenologs_env.out.env_path, // Needs both random and real data calculations to complete
                     compute_fdr_data.out.project_path,
                     compute_fdr_data.out.random_trials_sig,
                     compute_real_phenolog_data.out.real_sig,
                     taxon_id,
                     prd)

    // Compute ortholog to phenotype rankings
    compute_ortholog_rank_calcs(get_phenologs_env.out.env_path,
                                compute_fdr_info.out.project_path, 
                                taxon_id,
                                prd,
                                fdr,
                                kneighbs,
                                rank_metric)

    // Covert to similarity tables for alternate downstream processing
    convert_to_sim_tables(get_phenologs_env.out.env_path,
                          compute_ortholog_rank_calcs.out.project_path,
                          taxon_id,
                          prd,
                          fdr)

    // Leave one out cross validation
    if (params.xvalidate_calc) {
        leave_one_out_calculations(get_phenologs_env.out.env_path,
                                   compute_fdr_info.out.project_path,
                                   taxon_id,
                                   prd,
                                   fdr)
                                   } 

    // Leave one out cross validation assessment
    if (params.xvalidate_rank & params.xvalidate_calc) {
        leave_one_out_ortholog_rank_calcs(get_phenologs_env.out.env_path, 
                                          leave_one_out_calculations.out.project_path,
                                          taxon_id,
                                          prd,
                                          fdr,
                                          kneighbs,
                                          rank_metric)
                                          }
                                          
}