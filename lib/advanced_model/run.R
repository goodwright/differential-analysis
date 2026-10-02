#!/usr/bin/env Rscript
# Advanced mode: one fitted voom model, explicit numeric contrasts, optional CAMERA.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("Expected: options.json counts samples contrasts gene_sets universe")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "helpers.R"))
suppressPackageStartupMessages({library(jsonlite); library(edgeR); library(limma)})
opt <- fromJSON(args[1])
set.seed(as.integer(opt$seed))
for (key in c("min_samples", "min_set_size", "max_set_size")) {
    assert(length(opt[[key]]) == 1 && is.finite(opt[[key]]) && opt[[key]] == floor(opt[[key]]) && opt[[key]] >= 1,
           paste("Expected positive integer option:", key))
}
assert(opt$min_set_size >= 2 && opt$max_set_size >= opt$min_set_size, "Module sizes must satisfy 2 <= minimum <= maximum")
assert(length(opt$min_cpm) == 1 && is.finite(opt$min_cpm) && opt$min_cpm >= 0, "Minimum CPM must be finite and nonnegative")
assert(opt$set_adjust_method %in% c("BH", "holm", "bonferroni"), "Unsupported module adjustment method")
samples <- read_table(args[3])
original_sample_order <- samples$sample_id
# Match the existing Flow samplesheet checker's deterministic sample order.
samples <- samples[order(samples$sample_id), , drop = FALSE]
model <- validated_design(samples, opt$design_formula, opt$numeric_covariates)
design <- model$matrix
write_tsv(data.frame(sample_id = rownames(design), design, check.names = FALSE), "design_matrix.tsv")
contrast <- validated_contrasts(read_table(args[4]), colnames(design))
write_tsv(data.frame(coefficient = rownames(contrast), contrast, check.names = FALSE), "contrast_matrix.tsv")
raw <- read_table(args[2])
assert("gene_id" %in% names(raw), "Counts require a gene_id column")
assert(!anyDuplicated(names(raw)), "Duplicate count matrix columns")
assert(!anyNA(raw$gene_id) && all(nzchar(raw$gene_id)) && !anyDuplicated(raw$gene_id), "Gene IDs must be present and unique; aggregate intentionally upstream")
ids <- as.character(samples$sample_id)
assert(all(ids %in% names(raw)), "Counts are missing samples from the samplesheet")
assert(all(names(raw) %in% c("gene_id", "gene_name", ids)), "Unexpected count columns: supply exactly the selected samples, gene_id and optional gene_name")
counts <- as.matrix(raw[, ids, drop = FALSE]); storage.mode(counts) <- "numeric"
assert(all(is.finite(counts)) && all(counts >= 0), "Counts must be finite and nonnegative")
assert(all(colSums(counts) > 0), "Zero-depth library")
assert(opt$count_matrix_type %in% c("raw_counts", "length_scaled_counts"), "Unsupported count matrix type")
if (opt$count_matrix_type == "raw_counts") assert(all(abs(counts - round(counts)) < 1e-6), "Fractional counts require explicit length_scaled_counts mode; TPM/log expression is not supported")
rownames(counts) <- raw$gene_id
filter_samples <- rep(TRUE, nrow(samples))
if (nzchar(opt$filter_column)) {
    assert(opt$filter_column %in% names(samples), "Unknown expression-filter column")
    filter_samples <- samples[[opt$filter_column]] == opt$filter_level
}
assert(sum(filter_samples) >= opt$min_samples, "Too few samples for the prespecified expression filter")
keep <- rowSums(edgeR::cpm(counts[, filter_samples, drop = FALSE]) >= opt$min_cpm) >= opt$min_samples
id_column <- opt$gene_set_id_column
assert(id_column %in% names(raw), "Gene-set identifier column missing from counts")
gene_ids <- as.character(raw[[id_column]])
assert(!anyNA(gene_ids) && all(nzchar(gene_ids)) && !anyDuplicated(gene_ids), "Gene-set IDs must be nonempty and one-to-one; resolve aliases/duplicates upstream")
sets <- if (opt$use_gene_sets) read_gmt(args[5]) else list()
universe <- if (opt$use_universe) unique(trimws(readLines(args[6], warn = FALSE))) else gene_ids
assert(length(universe) > 0 && all(nzchar(universe)), "Empty gene universe")
if (opt$use_universe) assert(all(universe %in% gene_ids), "Universe contains IDs absent from the count matrix")
if (length(sets)) {
    absent_from_universe <- setdiff(intersect(unique(unlist(sets)), gene_ids), universe)
    assert(length(absent_from_universe) == 0, "Universe must include all observed module members; cannot selectively remove module genes")
}
write_tsv(data.frame(gene_id = raw$gene_id, set_id = gene_ids, passes_expression = keep,
                     in_universe = gene_ids %in% universe), "gene_filter_audit.tsv")
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
for (name in colnames(contrast)) {
    tab <- topTable(fit, coef = name, number = Inf, sort.by = "none", confint = TRUE, adjust.method = "BH")
    write_tsv(data.frame(gene_id = rownames(tab), set_id = retained_ids, tab), paste0(name, ".limma.results.tsv"))
}
# Correct across all requested module-by-contrast tests in this run; retain raw p-values.
camera_rows <- list()
for (name in colnames(contrast)) {
    if (!length(indexes)) next
    cam <- camera(v, index = indexes, design = design, contrast = contrast[, name],
                  inter.gene.cor = NA, allow.neg.cor = FALSE, sort = FALSE)
    cam$gene_set <- rownames(cam); cam$contrast <- name
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
write_json(list(options = opt, formula = model$formula, factor_levels = model$factor_levels,
    residual_df = model$residual_df, original_sample_order = original_sample_order, samples = ids, coefficients = colnames(design),
    genes_before_filter = nrow(counts), genes_after_expression_filter = sum(keep),
    genes_in_competitive_universe = nrow(v),
    module_correlation = "CAMERA estimates residual correlation per set/contrast, retaining finite residual degrees of freedom; negative estimates are truncated for variance inflation",
    normalization = "TMM on all expression-eligible genes before restricting the competitive universe",
    set_p_adjustment_family = "All tested gene-set-by-contrast combinations in this execution",
    gene_p_adjustment_family = "All retained genes, separately for each contrast"), "model_audit.json", pretty = TRUE, auto_unbox = TRUE)
saveRDS(list(voom = v, fit = fit, design = design, contrasts = contrast, options = opt), "fitted_model.rds")
capture.output(sessionInfo(), file = "R_sessionInfo.log")
writeLines(c('R_CUSTOM_MODEL:', paste0('    R: ', getRversion()),
    paste0('    limma: ', packageVersion('limma')), paste0('    edgeR: ', packageVersion('edgeR')),
    paste0('    statmod: ', packageVersion('statmod')), paste0('    jsonlite: ', packageVersion('jsonlite'))), "versions.yml")
