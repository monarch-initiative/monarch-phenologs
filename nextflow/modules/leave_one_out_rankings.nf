#!/usr/bin/env nextflow

process leave_one_out_ortholog_rank_calcs {
    tag 'leave_one_out_calculations'
    publishDir "./", mode: 'copy'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val taxon_id
    val prd
    val fdr
    val kneighbs
    val rank_metric

    output:
    path phenologs_data_dir, emit: project_path

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/phenologs_xvalidation_ortholog_to_phenotype_ranking.py -p ${phenologs_data_dir} \
                                                                                           -c ${task.cpus} \
                                                                                           -taxon_id ${taxon_id} \
                                                                                           -prd ${prd} \
                                                                                           -fdr ${fdr} \
                                                                                           -kneighbs ${kneighbs} \
                                                                                           -rank_metric ${rank_metric}
    """
}