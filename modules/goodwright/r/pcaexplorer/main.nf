process R_PCAEXPLORER {
    tag "$meta"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bioconductor-pcaexplorer=2.24.0" : null)
    container 'chrischeshire/pcaexplorer:latest'

    input:
    tuple val(meta), path(rds)
    val contrast
    val reference
    val treatment
    val blocking

    output:
    tuple val(meta), path("*.png")                , emit: png
    tuple val(meta), path("*.R_sessionInfo.log")  , emit: session_info
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    contrast_variable = contrast ?: "condition"
    reference_level = reference ?: "NO_REF"
    treatment_level = treatment ?: "NO_TREAT"
    blocking_variables = blocking ?: ""
    template 'pcaexplorer.R'

}
