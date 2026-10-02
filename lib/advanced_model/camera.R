#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "prepare.R"))
suppressPackageStartupMessages(library(limma))
assert(opt$use_gene_sets, "CAMERA requires fixed gene sets")
# Normalize using every expression-eligible gene, before restricting the competitive universe.
assert(sum(keep) >= 20, "Fewer than 20 expression-eligible genes")
y <- DGEList(counts[keep, , drop = FALSE]); y <- calcNormFactors(y, method = "TMM")
v <- voom(y, design, plot = FALSE)
write_tsv(data.frame(sample_id = ids, y$samples), "normalisation_factors.tsv")
retained_ids <- gene_ids[keep]
in_universe <- retained_ids %in% universe
assert(sum(in_universe) >= 20, "Fewer than 20 expression-eligible genes in the competitive universe")
v <- v[in_universe, ]
retained_ids <- retained_ids[in_universe]
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit, robust = TRUE)
write_tsv(data.frame(gene_id = rownames(v), set_id = retained_ids, v$E, check.names = FALSE), "voom_logcpm.tsv")
coverage <- list(); indexes <- list()
for (name in names(sets)) {
    members <- sets[[name]]; found <- match(members, retained_ids); found <- found[!is.na(found)]
    eligible <- length(found) >= opt$min_set_size && length(found) <= opt$max_set_size
    assert(length(found) < nrow(v), paste("Gene set has no competitive background:", name))
    coverage[[name]] <- data.frame(gene_set = name, requested = length(members), eligible = length(found),
        missing_or_filtered = paste(setdiff(members, retained_ids), collapse = ";"), tested = eligible)
    if (eligible) indexes[[name]] <- found
}
if (length(coverage)) write_tsv(do.call(rbind, coverage), "gene_set_coverage.tsv")
# Correct across all requested module-by-contrast tests in this run; retain raw p-values.
camera_rows <- list()
for (name in colnames(contrast)) {
    if (!length(indexes)) next
    cam <- camera(v, index = indexes, design = design, contrast = contrast[, name],
                  inter.gene.cor = NA, allow.neg.cor = FALSE, sort = FALSE)
    cam$gene_set <- rownames(cam); cam$contrast <- name
    cam$effect_estimator <- "limma-voom; not DESeq2"
    cam$FDR <- NULL
    cam$mean_set_log2_change <- vapply(indexes, function(ix) mean(fit$coefficients[ix, name]), 0.0)
    cam$mean_background_log2_change <- vapply(indexes, function(ix) mean(fit$coefficients[-ix, name]), 0.0)
    cam$excess_log2_change <- cam$mean_set_log2_change - cam$mean_background_log2_change
    cam$fraction_negative <- vapply(indexes, function(ix) mean(fit$coefficients[ix, name] < 0), 0.0)
    camera_rows[[name]] <- cam
}
if (length(camera_rows)) {
    result <- do.call(rbind, camera_rows)
    result$adjusted_p <- p.adjust(result$PValue, method = opt$set_adjust_method)
    write_tsv(result, "camera.results.tsv")
} else if (length(sets)) {
    write_tsv(data.frame(status = "No gene set passed the prespecified coverage thresholds; no module test performed"), "camera.not_tested.tsv")
}
write_json(list(engine = "limma-voom CAMERA", gene_results_source = "Separate DESeq2 process; not used as CAMERA input", options = opt, formula = model$formula, factor_levels = model$factor_levels,
    residual_df = model$residual_df, original_sample_order = original_sample_order, samples = ids, coefficients = colnames(design),
    genes_before_filter = nrow(counts), genes_after_expression_filter = sum(keep),
    genes_in_competitive_universe = nrow(v),
    module_correlation = "CAMERA estimates residual correlation per set/contrast, retaining finite residual degrees of freedom; negative estimates are truncated for variance inflation",
    normalization = "TMM on all expression-eligible genes before restricting the competitive universe",
    set_p_adjustment_family = "All tested gene-set-by-contrast combinations in this execution",
    gene_p_adjustment_family = "All retained genes, separately for each contrast"), "model_audit.json", pretty = TRUE, auto_unbox = TRUE)
saveRDS(list(voom = v, fit = fit, design = design, contrasts = contrast, options = opt), "fitted_model.rds")
capture.output(sessionInfo(), file = "R_sessionInfo.log")
writeLines(c('R_CAMERA_MODULE:', paste0('    R: ', getRversion()),
    paste0('    limma: ', packageVersion('limma')), paste0('    edgeR: ', packageVersion('edgeR')),
    paste0('    statmod: ', packageVersion('statmod')), paste0('    jsonlite: ', packageVersion('jsonlite'))), "versions.yml")
