process R_VOLCANO_PLOT {
    tag "${contrast}:${reference}_${treatment}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bioconductor-degreport=1.34.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-degreport:1.34.0--r42hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-degreport:1.34.0--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(deseq_results)
    val contrast
    val reference
    val treatment
    val blocking
    val fold_change
    val p_value

    output:
    tuple val(meta), path("*.pdf")               , emit: pdf
    tuple val(meta), path("*.R_sessionInfo.log") , emit: session_info
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    contrast_variable = contrast ?: "condition"
    reference_level = reference ?: "NO_REF"
    treatment_level = treatment ?: "NO_TREAT"
    blocking_variables = blocking ?: ""
    template 'volcano_plot.R'
}
