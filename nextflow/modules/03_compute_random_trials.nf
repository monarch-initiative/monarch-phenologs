#!/usr/bin/env nextflow

process compute_fdr_data {
    tag 'compute_fdr_data'
    publishDir path "./", mode: 'copy'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val taxon_id
    val prd
    val n_random_trials

    output:
    path "${phenologs_data_dir}", emit: project_path
    val "done", emit: random_trials_sig

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/03_compute_phenologs_randomized_trials.py -p ${phenologs_data_dir} \
                                                                              -n ${n_random_trials} \
                                                                              -c ${task.cpus} \
                                                                              -taxon_id ${taxon_id} \
                                                                              -prd ${prd}
    """
}
