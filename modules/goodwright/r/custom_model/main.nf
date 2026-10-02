process R_DESIGN_DESEQ2 {
    tag 'design_contrasts'
    label 'process_medium'

    // Release operators supply a built, immutable image digest. Local tests use their pinned R environment.
    container params.custom_model_container

    input:
    path counts, stageAs: 'counts.input'
    path samplesheet, stageAs: 'samplesheet.csv'
    path contrasts, stageAs: 'contrasts.input'
    path gene_sets, stageAs: 'gene_sets.gmt'
    path universe, stageAs: 'universe.txt'
    val settings_b64
    path model_code

    output:
    path '*.deseq2.results.tsv', emit: results
    path '*.tsv', emit: tables
    path 'model_audit.json', emit: audit
    path 'fitted_model.rds', emit: fit
    path 'R_sessionInfo.log', emit: session_info
    path 'versions.yml', emit: versions

    script:
    """
    printf '%s' '${settings_b64}' | base64 -d > model-options.json
    Rscript '${model_code}/run.R' model-options.json counts.input samplesheet.csv contrasts.input gene_sets.gmt universe.txt
    """
}

process R_CAMERA_MODULE {
    tag 'design_contrasts'
    label 'process_medium'

    // Release operators supply a built, immutable image digest. Local tests use their pinned R environment.
    container params.custom_model_container

    input:
    path counts, stageAs: 'counts.input'
    path samplesheet, stageAs: 'samplesheet.csv'
    path contrasts, stageAs: 'contrasts.input'
    path gene_sets, stageAs: 'gene_sets.gmt'
    path universe, stageAs: 'universe.txt'
    val settings_b64
    path model_code

    output:
    path 'camera.results.tsv', emit: camera, optional: true
    path 'camera.not_tested.tsv', emit: camera_not_tested, optional: true
    path 'gene_set_coverage.tsv', emit: coverage, optional: true
    path '*.tsv', emit: tables
    path 'model_audit.json', emit: audit
    path 'fitted_model.rds', emit: fit
    path 'R_sessionInfo.log', emit: session_info
    path 'versions.yml', emit: versions

    script:
    """
    printf '%s' '${settings_b64}' | base64 -d > model-options.json
    Rscript '${model_code}/camera.R' model-options.json counts.input samplesheet.csv contrasts.input gene_sets.gmt universe.txt
    """
}
