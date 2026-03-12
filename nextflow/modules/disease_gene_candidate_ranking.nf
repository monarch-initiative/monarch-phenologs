#!/usr/bin/env nextflow

process disease_gene_candidate_ranking {
    tag 'disease_gene_candidate_ranking'

    input:
    path phenologs_data_dir
    val taxon_id
    val prd
    val fdr
    val kneighbs
    val rank_metric

    output:
    path phenologs_data_dir, emit: project_path

    script:
    """
    python /app/python/phenologs_disease_gene_candidate_rankings.py -p ${phenologs_data_dir} \
                                                                     -taxon_id ${taxon_id} \
                                                                     -prd ${prd} \
                                                                     -fdr ${fdr} \
                                                                     -kneighbs ${kneighbs} \
                                                                     -rank_metric ${rank_metric}
    """
}
