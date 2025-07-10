#!/usr/bin/env nextflow

process convert_to_sim_tables {
    tag 'convert_to_sim_tables'
    publishDir "./", mode: 'copy'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val taxon_id
    val prd
    val fdr

    output:
    path phenologs_data_dir, emit: project_path

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/phenologs_to_oaksim.py -p ${phenologs_data_dir} \
                                                           -taxon_id ${taxon_id} \
                                                           -prd ${prd} \
                                                           -fdr ${fdr}
    """

}