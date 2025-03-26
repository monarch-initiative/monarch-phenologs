#!/usr/bin/env nextflow

process compute_ortholog_rank_calcs {
    tag 'compute_ortholog_rank_calcs'
    publishDir "./", mode: 'copy', overwrite: true

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val taxon_id
    val prd
    val fdr
    val kneighbs
    val rank_metric

    output:
    path phenologs_data_dir, emit: project_path
    path "monarch-phenologs", emit: project_env

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/06_phenologs_compute_ortholog_to_phenotype_rankings.py -p ${phenologs_data_dir} \
                                                                                           -taxon_id ${taxon_id} \
                                                                                           -prd ${prd} \
                                                                                           -fdr ${fdr} \
                                                                                           -kneighbs ${kneighbs} \
                                                                                           -rank_metric ${rank_metric}
    """

}
