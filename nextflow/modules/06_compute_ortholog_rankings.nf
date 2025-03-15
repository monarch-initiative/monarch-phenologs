#!/usr/bin/env nextflow

process compute_ortholog_rank_calcs {
    tag 'compute_ortholog_rank_calcs'
    publishDir "./", mode: 'copy'

    cpus 1
    memory '12GB'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val comp_sig_fdr
    val params

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_ranks

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/06_phenologs_compute_ortholog_to_phenotype_rankings.py -p ${phenologs_data_dir} \
                                                                                           -taxon_id ${params.taxon_id} \
                                                                                           -prd ${params.prd} \
                                                                                           -fdr ${params.fdr} \
                                                                                           -kneighbs ${params.kneighbs} \
                                                                                           -rank_metric ${params.rank_metric}
    """

}
