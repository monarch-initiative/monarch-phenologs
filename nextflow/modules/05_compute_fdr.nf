#!/usr/bin/env nextflow

process compute_fdr_info {
    tag 'compute_phenolog_data'
    publishDir path "./", mode: 'copy'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val random_trials_sig
    val real_sig
    val taxon_id
    val prd

    output:
    path "${phenologs_data_dir}", emit: project_path

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/05_compute_phenologs_fdr_data.py -p ${phenologs_data_dir} \
                                                                     -taxon_id ${taxon_id} \
                                                                     -prd ${prd}
    """
}
