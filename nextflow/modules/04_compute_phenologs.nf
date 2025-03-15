#!/usr/bin/env nextflow

process compute_real_phenolog_data {
    tag 'compute_phenolog_data'
    publishDir "./", mode: 'copy'
    
    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val params

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_real

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/04_compute_phenologs.py -p ${phenologs_data_dir} \
                                                            -c ${params.cpu_cores} \
                                                            -taxon_id ${params.taxon_id} \
                                                            -prd ${params.prd}
    """
}