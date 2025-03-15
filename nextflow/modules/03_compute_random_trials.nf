#!/usr/bin/env nextflow

process compute_fdr_data {
    tag 'compute_fdr_data'
    publishDir "./", mode: 'copy'

    cpus "${params.cpu_cores}"
    memory '64GB'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val params

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_trials

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/03_compute_phenologs_randomized_trials.py -p ${phenologs_data_dir} \
                                                                              -n ${params.n_random_trials} \
                                                                              -c ${params.cpu_cores} \
                                                                              -taxon_id ${params.taxon_id} \
                                                                              -prd ${params.prd}
    """
}
