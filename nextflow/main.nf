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

process compute_fdr_data {
    tag 'compute_fdr_data'
    publishDir "./", mode: 'copy'
    cpus 10
    memory '64GB'
    input:
    path phenologs_env_dir
    path phenologs_data_dir

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_trials

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/03_compute_phenologs_randomized_trials.py -p ${phenologs_data_dir} \
                                                                              -n ${params.n_random_trials} \
                                                                              -c ${params.cpu_cores} \
                                                                              -taxon_id ${params.taxon_id} \
                                                                              -prd ${params.prd}
    """
}

process compute_real_phenolog_data {
    tag 'compute_phenolog_data'
    publishDir "./", mode: 'copy'
    cpus 10
    memory '64GB'
    input:
    path phenologs_env_dir
    path phenologs_data_dir

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_real

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/04_compute_phenologs.py -p ${phenologs_data_dir} \
                                                            -c ${params.cpu_cores} \
                                                            -taxon_id ${params.taxon_id} \
                                                            -prd ${params.prd}
    """
}

process compute_fdr_info {
    tag 'compute_phenolog_data'
    publishDir "./", mode: 'copy'
    cpus 10
    memory '64GB'
    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val comp_sig_trials
    val comp_sig_real

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_fdr

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/05_compute_phenologs_fdr_data.py -p ${phenologs_data_dir} \
                                                                     -c ${params.cpu_cores} \
                                                                     -taxon_id ${params.taxon_id} \
                                                                     -prd ${params.prd}
    """
}

process compute_ortholog_rank_calcs {
    tag 'compute_ortholog_rank_calcs'
    publishDir "./", mode: 'copy'
    cpus 10
    memory '64GB'
    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val comp_sig_fdr

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_ranks

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/06_phenologs_compute_ortholog_to_phenotype_rankings.py -p ${phenologs_data_dir} \
                                                                                           -taxon_id ${params.taxon_id} \
                                                                                           -prd ${params.prd} \
                                                                                           -fdr ${params.fdr} \
                                                                                           -kneighbs ${params.kneighbs} \
                                                                                           -rank_metric ${params.rank_metric}
    """

}

process leave_one_out_calculations {
    tag 'leave_one_out_calculations'
    publishDir "./", mode: 'copy'
    cpus 10
    memory '64GB'
    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val comp_sig_fdr

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_xval

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/phenologs_leave_one_out_xvalidation.py -p ${phenologs_data_dir} \
                                                                           -c ${params.cpu_cores} \
                                                                           -taxon_id ${params.taxon_id} \
                                                                           -prd ${params.prd} \
                                                                           -fdr ${params.fdr}
    """
}

process leave_one_out_ortholog_rank_calcs {
    tag 'leave_one_out_calculations'
    publishDir "./", mode: 'copy'
    cpus 10
    memory '64GB'
    input:
    path phenologs_env_dir
    path phenologs_data_dir
    val comp_sig_xval

    output:
    path "phenologs-from-kg", emit: project_path
    val "done", emit: comp_sig_xval_rank

    script:
    """
    source ${phenologs_env_dir}/.venv/bin/activate
    python monarch-phenologs/python/phenologs_xvalidation_ortholog_to_phenotype_ranking.py -p ${phenologs_data_dir} \
                                                                                           -c ${params.cpu_cores} \
                                                                                           -taxon_id ${params.taxon_id} \
                                                                                           -prd ${params.prd} \
                                                                                           -fdr ${params.fdr} \
                                                                                           -kneighbs ${params.kneighbs} \
                                                                                           -rank_metric ${params.rank_metric}
    """
}

// Phenologs calculation parameters 
params.n_random_trials = 2
params.cpu_cores = 10
params.taxon_id = 9606
params.prd = "disease"

// Phenologs ortholog to phenotype ranking params
params.fdr = .95
params.kneighbs = 10
params.rank_metric = 'nb'

// Xvalidation options
params.xvalidate_calc = false
params.xvalidate_rank = false


workflow {

    // Get the monarch-phenologs environment
    get_phenologs_env()
    get_phenologs_data(get_phenologs_env.out)

    // Compute random trials
    compute_fdr_data(get_phenologs_env.out, 
                     get_phenologs_data.out)
    
    // Compute real phenolog data
    compute_real_phenolog_data(get_phenologs_env.out, 
                               get_phenologs_data.out)
    
    // Compute FDR info for significance cut
    compute_fdr_info(get_phenologs_env.out, 
                     get_phenologs_data.out,
                     compute_fdr_data.out.comp_sig_trials,          // Completion signal for random trials
                     compute_real_phenolog_data.out.comp_sig_real)  // Completion signal for real phenolog calcs

    // Compute ortholog to phenotype rankings
    compute_ortholog_rank_calcs(get_phenologs_env.out, 
                                get_phenologs_data.out,
                                compute_fdr_info.out.comp_sig_fdr)  // Completion signal for FDR info 
    
    // Leave one out cross validation
    if (params.xvalidate_calc) {
        leave_one_out_calculations(get_phenologs_env.out, 
                                   get_phenologs_data.out,
                                   compute_fdr_info.out.comp_sig_fdr)
                                   }  // Completion signal for FDR info
    
    // Leave one out cross validation assessment
    if (params.xvalidate_rank) {
        leave_one_out_ortholog_rank_calcs(get_phenologs_env.out, 
                                        get_phenologs_data.out,
                                        leave_one_out_calculations.out.comp_sig_xval)  // Completion signal for xval calcs
    }

}