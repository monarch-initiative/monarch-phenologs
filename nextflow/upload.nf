#!/usr/bin/env nextflow

include { upload_results } from './modules/upload_results.nf'

params.base_dir = false
params.release_dir = false


workflow {
    // Call our upload module
    upload_results(params.base_dir, params.release_dir)
}
