#!/usr/bin/env nextflow

process get_phenologs_data {
    tag "get_phenologs_data"
    publishDir "./", mode: 'copy'

    input:
    path phenologs_env_dir

    output:
    path "phenologs-from-kg"

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/01_download_and_setup.py -p ./phenologs-from-kg
    python monarch-phenologs/python/02_kg_data_extraction.py -p ./phenologs-from-kg
    cd phenologs-from-kg
    """
}