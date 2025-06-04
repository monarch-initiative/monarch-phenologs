#!/usr/bin/env nextflow

process merge_outputs {
    tag "merge_outputs"

    input:
    path base_dir
    path random_results
    path real_results

    output:
    // Note, that we are using base_dir as a variable name rather than overwriting the input base_dir
    path base_dir, emit: project_path

    script:
    """
    """
}