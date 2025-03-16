#!/usr/bin/env nextflow

process compute_fdr_data {
    tag 'compute_fdr_data'
    publishDir "./", mode: 'copy'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val params

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: random_sig

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/03_compute_phenologs_randomized_trials.py -p ${phenologs_data_dir} \
                                                                              -n ${params.n_random_trials} \
                                                                              -c ${task.cpus} \
                                                                              -taxon_id ${params.taxon_id} \
                                                                              -prd ${params.prd}
    """
}
