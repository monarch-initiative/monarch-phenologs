#!/usr/bin/env nextflow

process get_phenologs_data {
    tag "get_phenologs_data"

    input:
    val graph_version
    val phenio_version

    output:
    path "phenologs-from-kg", emit: project_path

    script:
    """
    python /app/python/01_download_and_setup.py -p ./phenologs-from-kg \
                                                 -gv ${graph_version} \
                                                 -pv ${phenio_version}
    python /app/python/02_kg_data_extraction.py -p ./phenologs-from-kg
    """
}
