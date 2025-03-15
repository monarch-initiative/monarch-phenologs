#!/usr/bin/env nextflow

process get_phenologs_env {
    tag "get_phenologs_env"
    publishDir "./", mode: 'copy'

    output:
    path "monarch-phenologs"

    script:
    """
    git clone -b monarch_phenologs_optimized_fdr_calcs git@github.com:monarch-initiative/monarch-phenologs.git
    cd monarch-phenologs
    poetry config virtualenvs.in-project true
    poetry install
    """
}