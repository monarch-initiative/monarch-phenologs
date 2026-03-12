#!/usr/bin/env nextflow

process compute_fdr_data {
    tag 'compute_fdr_data'

    input:
    path phenologs_data_dir
    val taxon_id
    val prd
    val n_random_trials

    output:
    path "${phenologs_data_dir}/random_trials/*", emit: data_path

    script:
    """
    python /app/python/03_compute_phenologs_randomized_trials.py -p ${phenologs_data_dir} \
                                                                  -n ${n_random_trials} \
                                                                  -c ${task.cpus} \
                                                                  -taxon_id ${taxon_id} \
                                                                  -prd ${prd}
    """
}
