#!/usr/bin/env nextflow

process leave_one_out_calculations {
    tag 'leave_one_out_calculations'
    publishDir "./", mode: 'copy'

    cpus "${params.cpu_cores}"
    memory '64GB'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val comp_sig_fdr
    val params

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_xval

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/phenologs_leave_one_out_xvalidation.py -p ${phenologs_data_dir} \
                                                                           -c ${params.cpu_cores} \
                                                                           -taxon_id ${params.taxon_id} \
                                                                           -prd ${params.prd} \
                                                                           -fdr ${params.fdr}
    """
}
