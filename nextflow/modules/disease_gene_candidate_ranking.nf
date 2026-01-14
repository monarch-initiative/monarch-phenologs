#!/usr/bin/env nextflow

process disease_gene_candidate_ranking {
    tag 'disease_gene_candidate_ranking'

    input:
    path phenologs_data_dir
    val fdr
    val rank_type

    output:
    path phenologs_data_dir, emit: project_path

    script:
    """
    python /app/python/phenologs_disease_gene_candidate_rankings.py -p ${phenologs_data_dir} \
                                                                     -fdr ${fdr} \
                                                                     -rank_type ${rank_type}
    """
}
