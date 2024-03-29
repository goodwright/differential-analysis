process R_PCAEXPLORER {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bioconductor-pcaexplorer=2.24.0" : null)
    container 'chrischeshire/pcaexplorer:latest'

    input:
    tuple val(meta), path(rds)
    val contrast
    val blocking

    output:
    tuple val(meta), path("*.pdf")             , emit: pdf
    tuple val(meta), path("R_sessionInfo.log") , emit: session_info
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    contrast_variable = contrast ?: "condition"
    blocking_variables = blocking ?: ""
    template 'pcaexplorer.R'

}
