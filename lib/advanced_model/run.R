#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "prepare.R"))
suppressPackageStartupMessages(library(DESeq2))
assert(opt$deseq_fit_type %in% c("parametric", "local", "mean"), "Unsupported DESeq2 dispersion fit")
assert(opt$deseq_sf_type %in% c("ratio", "poscounts"), "Unsupported DESeq2 size-factor method")
# DESeq2's tximport constructor also rounds count estimates. This conversion is explicit and audited.
model_counts <- counts[keep, , drop = FALSE]
rounding <- list(input_type = opt$count_matrix_type, fractional_entries = sum(model_counts != round(model_counts)),
                 maximum_absolute_change = max(abs(model_counts - round(model_counts))),
                 rule = "Nearest integer using R round; only explicitly declared lengthScaledTPM estimates may be fractional")
assert(all(model_counts <= .Machine$integer.max), "Counts exceed DESeq2 integer range")
model_counts <- round(model_counts)
storage.mode(model_counts) <- "integer"
assert(all(rowSums(model_counts) > 0), "Expression-eligible gene became all zero after count conversion")
metadata <- model$samples; rownames(metadata) <- ids
# Pass the already validated matrix so numeric contrasts use exactly its saved column order.
dds <- DESeqDataSetFromMatrix(countData = model_counts, colData = metadata, design = design)
dds <- DESeq(dds, test = "Wald", fitType = opt$deseq_fit_type, sfType = opt$deseq_sf_type,
             betaPrior = FALSE, minReplicatesForReplace = Inf, quiet = TRUE)
assert(all(mcols(dds)$betaConv), "DESeq2 coefficients did not converge; no contrasts were reported")
assert(length(resultsNames(dds)) == ncol(design), "DESeq2 coefficient count changed")
write_tsv(data.frame(design_coefficient = colnames(design), deseq_coefficient = resultsNames(dds)), "coefficient_mapping.tsv")
write_tsv(data.frame(sample_id = ids, size_factor = sizeFactors(dds)), "normalisation_factors.tsv")
write_tsv(data.frame(gene_id = rownames(dds), counts(dds, normalized = TRUE), check.names = FALSE), "normalised_counts.tsv")
for (name in colnames(contrast)) {
    result <- results(dds, contrast = unname(contrast[, name]), independentFiltering = FALSE,
                      pAdjustMethod = "BH", lfcThreshold = 0, altHypothesis = "greaterAbs")
    tab <- as.data.frame(result)
    tab$CI95_low <- tab$log2FoldChange - qnorm(.975) * tab$lfcSE
    tab$CI95_high <- tab$log2FoldChange + qnorm(.975) * tab$lfcSE
    tab$test_status <- ifelse(is.na(tab$pvalue), "not_tested_check_Cooks_distance", "tested")
    write_tsv(data.frame(gene_id = rownames(tab), set_id = gene_ids[keep], tab), paste0(name, ".deseq2.results.tsv"))
}
write_tsv(data.frame(gene_id = rownames(dds), as.data.frame(mcols(dds)), check.names = FALSE), "gene_diagnostics.tsv")
write_tsv(data.frame(gene_id = rownames(dds), assays(dds)[["cooks"]], check.names = FALSE), "cooks_distances.tsv")
write_json(list(engine = "DESeq2", options = opt, formula = model$formula, factor_levels = model$factor_levels,
    residual_df = model$residual_df, original_sample_order = original_sample_order, samples = ids,
    coefficients = colnames(design), deseq_coefficients = resultsNames(dds), rounding = rounding,
    genes_before_filter = nrow(counts), genes_after_expression_filter = sum(keep),
    normalization = paste("DESeq2", opt$deseq_sf_type, "size factors on all expression-eligible genes; module universe does not restrict gene tests"),
    requested_dispersion_fit = opt$deseq_fit_type, actual_dispersion_fit = attr(dispersionFunction(dds), "fitType"),
    test = "Wald numeric contrast; unshrunk log2 fold changes; 95% normal-approximation intervals",
    independent_filtering = FALSE, beta_prior = FALSE, automatic_outlier_replacement = FALSE,
    cooks_filter = "DESeq2 default Cook's cutoff retained; excluded p-values remain NA with diagnostics",
    gene_p_adjustment_family = "BH across nonmissing gene p-values among all expression-eligible genes, separately for each contrast",
    module_testing = "Separate optional limma-voom CAMERA process; not a DESeq2 module p-value"), "model_audit.json", pretty = TRUE, auto_unbox = TRUE)
saveRDS(list(dds = dds, design = design, contrasts = contrast, options = opt), "fitted_model.rds")
capture.output(sessionInfo(), file = "R_sessionInfo.log")
writeLines(c('R_DESIGN_DESEQ2:', paste0('    R: ', getRversion()),
    paste0('    DESeq2: ', packageVersion('DESeq2')), paste0('    edgeR: ', packageVersion('edgeR')),
    paste0('    jsonlite: ', packageVersion('jsonlite'))), "versions.yml")
