#!/usr/bin/env nextflow

process leave_one_out_calculations {
    tag 'leave_one_out_calculations'

    input:
    path phenologs_data_dir
    val taxon_id
    val prd
    val fdr

    output:
    path phenologs_data_dir, emit: project_path

    script:
    """
    python /app/python/phenologs_leave_one_out_xvalidation.py -p ${phenologs_data_dir} \
                                                               -c ${task.cpus} \
                                                               -taxon_id ${taxon_id} \
                                                               -prd ${prd} \
                                                               -fdr ${fdr}
    """
}
