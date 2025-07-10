#!/usr/bin/env nextflow

process disease_gene_candidate_ranking {
    tag 'disease_gene_candidate_ranking'

    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val fdr

    output:
    path phenologs_data_dir, emit: project_path

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/phenologs_disease_gene_candidate_rankings.py -p ${phenologs_data_dir} \
                                                                                 -fdr ${fdr}
    """
}
