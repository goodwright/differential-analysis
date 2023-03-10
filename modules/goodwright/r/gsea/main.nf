process R_GSEA {
    tag "${meta.id}"
    label 'process_medium'

    container 'chrischeshire/gsea:latest'

    input:
    tuple val(meta), path(deseq2_results)
    val organism

    output:
    tuple val(meta), path("*.genekitr_gsea_result.xlsx"), emit: ora_data
    tuple val(meta), path("*.genekitr_ora_result.xlsx") , emit: gsea_data
    tuple val(meta), path("*.png")                      , emit: plots
    tuple val(meta), path("*.R_sessionInfo.log")        , emit: session_info
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'gsea.R'
}
