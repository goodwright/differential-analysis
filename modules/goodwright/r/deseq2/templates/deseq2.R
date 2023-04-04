#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse
parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- unlist(lapply(args_list, function(y) strsplit(y, ' +')), recursive = FALSE)

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Passed to read.delim()
#' @param row.names Passed to read.delim()
#'
#' @return output Data frame
read_delim_flexible <- function(file, header = TRUE, row.names = NULL){

    ext <- tolower(tail(strsplit(basename(file), split = "\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = FALSE
    )
}

#' Round numeric dataframe columns to fixed decimal places by applying
#' formatting and converting back to numerics
#'
#' @param dataframe A data frame
#' @param columns Which columns to round (assumes all of them by default)
#' @param digits How many decimal places to round to?
#'
#' @return output Data frame
round_dataframe_columns <- function(df, columns = NULL, digits = 8){
    if (is.null(columns)){
        columns <- colnames(df)
    }

    df[,columns] <- format(data.frame(df[, columns]), nsmall = digits)

    # Convert columns back to numeric

    for (c in columns) {
        df[[c]][grep("^ *NA$", df[[c]])] <- NA
        df[[c]] <- as.numeric(df[[c]])
    }
    df
}

#####################################################
#####################################################
## PARSE PARAMETERS FROM NEXTFLOW AND COMMAND LINE ##
#####################################################
#####################################################

# Set defaults and classes
opt <- list(
    cores = 1,                                    # Number of cores to use
    count_file = "!{counts}",                     # The input count matrix file
    sample_file = '!{samplesheet}',               # The experimental design file
    contrast_variable = "!{contrast_variable}",   # The design column that will be used for comparision
    reference_level = "!{reference_level}",       # The reference level name for the contrast
    treatment_level = "!{treatment_level}",       # The treatment level name for the contrast
    blocking_variables = "!{blocking_variables}", # The design columns that will be used for blocking factors
    gene_id_col = "gene_id",                      # The fault id column in the count matrix
    sample_id_col = "sample_id",                  # The sample id column sample sheet

    # DESeq params, these are the defaults - see https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf
    test = "Wald",
    fit_type = "parametric",
    min_replicates_for_replace = 7,
    use_t = FALSE,
    sf_type = 'ratio',

    # results params, these are the defaults - see https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf
    lfc_threshold = 0,
    alt_hypothesis = 'greaterAbs',
    independent_filtering = TRUE,
    p_adjust_method = 'BH',
    alpha = 0.1,
    minmu = 0.5,

    shrink_lfc = TRUE,
    lfcshrink_type = 'ashr'

    # vs_method = 'vst', # 'rlog', 'vst', or 'rlog,vst'
    # vs_blind = TRUE,
    # vst_nsub = 1000
)
opt_types <- lapply(opt, class)

# Parse command line args
cl_args <- commandArgs(trailingOnly=TRUE)
cl_keys <- grep("^--", cl_args, value = TRUE)
cl_opt <- list()
for (key in cl_keys) {
    key_index <- which(cl_args == key)
    value <- cl_args[key_index + 1]
    cl_opt[[sub("^--", "", key)]] <- value
}

# Override defaults with command line args
opt <- modifyList(opt, cl_opt)

# Apply parameter overrides
args_opt <- parse_args("!{task.ext.args}")
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else{
        # Preserve classes from defaults where possible
        if (! is.null(opt[[ao]])){
            args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
        }
        opt[[ao]] <- args_opt[[ao]]
    }
}

if (is.na(opt$blocking_variables) || opt$blocking_variables== '') {
    opt$blocking_variables <- NULL
}

# Check if required parameters have been provided
required_opts <- c('contrast_variable', 'reference_level', 'treatment_level')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('count_file', 'sample_file')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

print(opt)

################################################
################################################
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
    read_delim_flexible(
        file = opt$count_file,
        header = TRUE,
        row.names = opt$gene_id_col
    )
sample.sheet <- read_delim_flexible(file = opt$sample_file)

# Prepare the sample sheet
rownames(sample.sheet) <- sample.sheet[[opt$sample_id_col]]

# head(count.table)

# Prepare the count table
# Save any non-count data, will gene metadata etc we might need later
# Round the count table values
noncount.table <- count.table[, !colnames(count.table) %in% rownames(sample.sheet), drop = FALSE]
count.table <- count.table[, rownames(sample.sheet)]
count.table <- round(count.table)

# ################################################
# ################################################
# ## CHECK MODEL SPECIFICATION                  ##
# ################################################
# ################################################

blocking.vars <- c()

if (!opt$contrast_variable %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable \"',
        opt$contrast_variable,
        '\" not in sample sheet'
        )
    )
} else if (any(!c(opt$reference_level, opt$treatment_level) %in% sample.sheet[[opt$contrast_variable]])) {
    stop(
        paste(
        'Please choose reference and treatment levels that are present in the',
        opt$contrast_variable,
        'column of the sample sheet'
        )
    )
} else if (!is.null(opt$blocking_variables)) {
    blocking.vars = unlist(strsplit(opt$blocking_variables, split = ';'))
    if (!all(blocking.vars %in% colnames(sample.sheet))) {
        missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
        stop(
            paste(
                'Blocking variables', missing_block,
                'do not correspond to sample sheet columns.'
            )
        )
    }
}

# Make sure all the appropriate variables are factors
for (v in c(blocking.vars, opt$contrast_variable)) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
}

# Now specify the model. Use cell-means style so we can be explicit with the
# contrasts
model <- '~ 0 +'

# ~ sex + condition
# ~ sex + condition + batch
# ~ sex + condition + batch + genotype + genotype:condition
# Variable of interest goes last, see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs
if (!is.null(opt$blocking_variables)) {
    model <- paste(model, paste(blocking.vars, collapse = '+'))
    model <- paste(model, opt$contrast_variable, sep = ' + ')
} else {
    model <- paste(model, opt$contrast_variable, sep = ' ')
}

model <- as.formula(model)
print(model)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(DESeq2)
library(BiocParallel)

################################################
################################################
## Run DESeq2 processes                       ##
################################################
################################################

dds <- DESeqDataSetFromMatrix(
    countData = count.table,
    colData = sample.sheet,
    design = model
)

dds <- DESeq(
    dds,
    test = opt$test,
    fitType = opt$fit_type,
    minReplicatesForReplace = opt$min_replicates_for_replace,
    useT = opt$use_t,
    sfType = opt$sf_type,
    parallel=TRUE, BPPARAM=MulticoreParam(opt$cores)
)

comp.results <-
    results(
        dds,
        lfcThreshold = opt$lfc_threshold,
        altHypothesis = opt$alt_hypothesis,
        independentFiltering = opt$independent_filtering,
        alpha = opt$alpha,
        pAdjustMethod = opt$p_adjust_method,
        minmu = opt$minmu,
        contrast = c(
            opt$contrast_variable,
            c(opt$treatment_level, opt$reference_level)
        )
    )

# Adds shrunken log2 fold changes (LFC) and SE to a results table from DESeq run without LFC shrinkage.
if (opt$shrink_lfc){
    comp.results <- lfcShrink(dds,
        type = opt$lfcshrink_type,
        contrast = c(
            opt$contrast_variable,
            c(opt$treatment_level, opt$reference_level)
        )
    )
}

################################################
################################################
## Generate outputs                           ##
################################################
################################################

prefix_part_names <- c('contrast_variable', 'reference_level', 'treatment_level', 'blocking_variables')
prefix_parts <- unlist(lapply(prefix_part_names, function(x) gsub("[^[:alnum:]]", "_", opt[[x]])))
output_prefix <- paste(prefix_parts[prefix_parts != ''], collapse = '-')

contrast.name <-
    paste(opt$treatment_level, opt$reference_level, sep = "_vs_")
cat("Saving results for ", contrast.name, " ...\n", sep = "")

# Differential expression table- note very limited rounding for consistency of
# results
write.table(
    data.frame(
        gene_id = rownames(comp.results),
        round_dataframe_columns(data.frame(comp.results))
    ),
    file = paste(output_prefix, 'deseq2.results.tsv', sep = '.'),
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# R objects for other processes to use
saveRDS(dds, file = paste(output_prefix, 'dds.rld.rds', sep = '.'))

# Size factors
sf_df = data.frame(sample = names(sizeFactors(dds)), data.frame(sizeFactors(dds)))
colnames(sf_df) <- c('sample', 'sizeFactor')
write.table(
    sf_df,
    file = 'deseq2.sizefactors.tsv',
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

# Write specified matrices
write.table(
    data.frame(gene_id=rownames(counts(dds)), counts(dds, normalized = TRUE)),
    file = 'normalised_counts.tsv',
    col.names = TRUE,
    row.names = FALSE,
    sep = '\t',
    quote = FALSE
)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
deseq2.version <- as.character(packageVersion('DESeq2'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    bioconductor-deseq2:', deseq2.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
