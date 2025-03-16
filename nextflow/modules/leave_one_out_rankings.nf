#!/usr/bin/env nextflow

process leave_one_out_ortholog_rank_calcs {
    tag 'leave_one_out_calculations'
    publishDir "./", mode: 'copy'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val params

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: xval_rank_sig

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/phenologs_xvalidation_ortholog_to_phenotype_ranking.py -p ${phenologs_data_dir} \
                                                                                           -c ${task.cpus} \
                                                                                           -taxon_id ${params.taxon_id} \
                                                                                           -prd ${params.prd} \
                                                                                           -fdr ${params.fdr} \
                                                                                           -kneighbs ${params.kneighbs} \
                                                                                           -rank_metric ${params.rank_metric}
    """
}