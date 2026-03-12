#!/usr/bin/env nextflow

process upload_results {
    tag 'upload_results'

    input:
    val base_dir
    val release_dir

    script:
    """
    ###################################
    ### Compress and upload results ###

    ### Change to our results directory
    cd ${base_dir}
    cd ${release_dir}

    ### Three files and a directory we want to upload for disease results
    DRES1="phenologs-from-kg/phenologs_results/Homo_sapiens_disease_results/Homo-sapiens_pooled_phenologs_fdr0.95.tsv"
    DRES2="phenologs-from-kg/phenologs_results/Homo_sapiens_disease_results/Homo-sapiens_xspecies_phenolog_counts.pdf"
    DRES3="phenologs-from-kg/phenologs_results/Homo_sapiens_disease_results/Homo-sapiens_xspecies_random_vs_real.pdf"
    DRES4="phenologs-from-kg/phenologs_results/Homo_sapiens_disease_results/Homo_sapiens_disease_gene_candidates_fdr0.95"

    ### Three files we want to upload for phenotype results
    PRES1="phenologs-from-kg/phenologs_results/Homo_sapiens_phenotype_results/Homo-sapiens_pooled_phenologs_fdr0.95.tsv"
    PRES2="phenologs-from-kg/phenologs_results/Homo_sapiens_phenotype_results/Homo-sapiens_xspecies_phenolog_counts.pdf"
    PRES3="phenologs-from-kg/phenologs_results/Homo_sapiens_phenotype_results/Homo-sapiens_xspecies_random_vs_real.pdf"

    ### We will upload base information we calculate as well for full transparency
    SPECIES_DATA="phenologs-from-kg/species_data"

    ### Google cloud bucket info we want to upload to
    UPLOAD_BASE="gs://data-public-monarchinitiative/phenologs/alpha"
    UPLOAD_BUCKET="\${UPLOAD_BASE}/${release_dir}"
    DISEASE_BUCKET="\${UPLOAD_BUCKET}/disease_results/"
    PHENOTYPE_BUCKET="\${UPLOAD_BUCKET}/phenotype_results/"

    ### Compress data where applicable
    gzip -c "\$DRES1" > "\${DRES1}.gz"
    gzip -c "\$PRES1" > "\${PRES1}.gz"
    tar -czf "\${DRES4}.tar.gz" "\$DRES4"
    tar -czf "\${SPECIES_DATA}.tar.gz" "\$SPECIES_DATA"

    ### Upload disease results
    gcloud storage cp "\${DRES1}.gz" "\$DISEASE_BUCKET"
    gcloud storage cp "\$DRES2" "\$DISEASE_BUCKET"
    gcloud storage cp "\$DRES3" "\$DISEASE_BUCKET"
    gcloud storage cp "\${DRES4}.tar.gz" "\$DISEASE_BUCKET"

    ### Upload phenotype results
    gcloud storage cp "\${PRES1}.gz" "\$PHENOTYPE_BUCKET"
    gcloud storage cp "\$PRES2" "\$PHENOTYPE_BUCKET"
    gcloud storage cp "\$PRES3" "\$PHENOTYPE_BUCKET"

    ### Upload input data to algorithm
    gcloud storage cp "\${SPECIES_DATA}.tar.gz" "\$UPLOAD_BUCKET"

    #####################
    ### Clean up step ###
    ###rm -rf work
    ###nextflow clean -f
    """
}
