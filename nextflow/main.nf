#!/usr/bin/env nextflow

include { get_phenologs_data } from './modules/02_get_data.nf'
include { compute_fdr_data } from './modules/03_compute_random_trials.nf'
include { compute_real_phenolog_data } from './modules/04_compute_phenologs.nf'
include { compute_fdr_info } from './modules/05_compute_fdr.nf'
include { compute_ortholog_rank_calcs } from './modules/06_compute_ortholog_rankings.nf'
include { convert_to_sim_tables } from './modules/to_similarity_tables.nf'
include { leave_one_out_calculations } from './modules/leave_one_out_calculations.nf'
include { leave_one_out_ortholog_rank_calcs } from './modules/leave_one_out_rankings.nf'
include { merge_outputs } from './modules/merge_outputs.nf'
include { disease_gene_candidate_ranking } from './modules/disease_gene_candidate_ranking.nf'


// Data download parameters
params.graph_version = "latest"   // Monarch KG version (e.g., "latest" or "2025-06-19")
params.phenio_version = "latest"  // Phenio relation graph version

// Phenologs calculation parameters
params.n_random_trials = 1
params.taxon_id = 9606
params.prd = "disease" // "phenotype" or "disease"

// Phenologs ortholog to phenotype ranking params
params.fdr = .95
params.kneighbs = 10
params.rank_metric = 'nb'
params.rank_type = 'protein_family'  // "protein_family" or "gene"

// Xvalidation options
params.xvalidate_calc = false
params.xvalidate_rank = false


workflow {

    // Set up the parameters for the workflow
    Channel.value(params.graph_version).set{ graph_version }
    Channel.value(params.phenio_version).set{ phenio_version }
    Channel.value(params.taxon_id).set{ taxon_id }
    Channel.value(params.prd).set{ prd }
    Channel.value(params.n_random_trials).set{ n_random_trials }
    Channel.value(params.fdr).set{ fdr }
    Channel.value(params.kneighbs).set{ kneighbs }
    Channel.value(params.rank_metric).set{ rank_metric }
    Channel.value(params.rank_type).set{ rank_type }

    // Get phenologs data (container provides the Python environment)
    get_phenologs_data(graph_version, phenio_version)

    // Run next two steps in parallel
    // Compute random trials
    compute_fdr_data(get_phenologs_data.out.project_path,
                     taxon_id,
                     prd,
                     n_random_trials)

    // Compute real phenolog data
    compute_real_phenolog_data(get_phenologs_data.out.project_path,
                               taxon_id,
                               prd)

    // Merge outputs from random trials and real phenolog data
    merge_outputs(get_phenologs_data.out.project_path,
                  compute_fdr_data.out.data_path,
                  compute_real_phenolog_data.out.data_path)


    // Compute FDR info for significance cut
    compute_fdr_info(merge_outputs.out.project_path,
                     taxon_id,
                     prd)

    // Compute ortholog to phenotype rankings
    compute_ortholog_rank_calcs(compute_fdr_info.out.project_path,
                                taxon_id,
                                prd,
                                fdr,
                                kneighbs,
                                rank_metric,
                                rank_type)

    // Covert to similarity tables for alternate downstream processing
    convert_to_sim_tables(compute_ortholog_rank_calcs.out.project_path,
                          taxon_id,
                          prd,
                          fdr)

    // Compute disease gene candidate rankings (only supports rank_type="gene")
    if (params.taxon_id == 9606 && params.prd == "disease") {
        disease_gene_candidate_ranking(convert_to_sim_tables.out.project_path,
                                       fdr,
                                       rank_type)
                                       }


    // Leave one out cross validation
    if (params.xvalidate_calc) {
        leave_one_out_calculations(compute_fdr_info.out.project_path,
                                   taxon_id,
                                   prd,
                                   fdr)
                                   }

    // Leave one out cross validation assessment
    if (params.xvalidate_rank & params.xvalidate_calc) {
        leave_one_out_ortholog_rank_calcs(leave_one_out_calculations.out.project_path,
                                          taxon_id,
                                          prd,
                                          fdr,
                                          kneighbs,
                                          rank_metric,
                                          rank_type)
                                          }
}
