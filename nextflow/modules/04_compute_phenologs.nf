#!/usr/bin/env nextflow

process compute_real_phenolog_data {
    tag 'compute_real_phenolog_data'
    
    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val taxon_id
    val prd

    output:
    path "${phenologs_data_dir}/phenologs_results/*", emit: data_path

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/04_compute_phenologs.py -p ${phenologs_data_dir} \
                                                            -c ${task.cpus} \
                                                            -taxon_id ${taxon_id} \
                                                            -prd ${prd}
    """
}