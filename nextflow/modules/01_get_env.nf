#!/usr/bin/env nextflow

process get_phenologs_env {
    tag "get_phenologs_env"

    output:
    path "monarch-phenologs", emit: env_path

    script:
    """
    git clone -b monarch_phenologs_optimized_fdr_calcs git@github.com:monarch-initiative/monarch-phenologs.git
    cd monarch-phenologs
    poetry config virtualenvs.in-project true
    poetry install
    """
}