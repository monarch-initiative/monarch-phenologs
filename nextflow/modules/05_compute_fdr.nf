#!/usr/bin/env nextflow

process compute_fdr_info {
    tag 'compute_phenolog_data'
    publishDir "./", mode: 'copy'

    input:
    path phenologs_env_dir
    path phenologs_data_dir // Both random and real data are needed and point to same directory (i think)
    val random_sig
    val real_sig
    val params

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: fdr_sig

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/05_compute_phenologs_fdr_data.py -p ${phenologs_data_dir} \
                                                                     -taxon_id ${params.taxon_id} \
                                                                     -prd ${params.prd}
    """
}
