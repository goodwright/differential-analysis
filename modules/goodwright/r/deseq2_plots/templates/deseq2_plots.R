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

#####################################################
#####################################################
## PARSE PARAMETERS FROM NEXTFLOW AND COMMAND LINE ##
#####################################################
#####################################################

# Set defaults and classes
opt <- list(
    deseq_rds = "!{rds}",
    contrast_variable = "!{contrast_variable}",
    reference_level = "!{reference_level}",
    treatment_level = "!{treatment_level}",
    blocking_variables = "!{blocking_variables}",

    plot_width = 1800,
    plot_height = 1200,
    plot_res = 300,
    plot_point_size = 4
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
required_opts <- c('deseq_rds', 'contrast_variable', 'reference_level', 'treatment_level')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('deseq_rds')){
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
## Load libraries                             ##
################################################
################################################

library(DESeq2)

################################################
################################################
## Load data and generate plots               ##
################################################
################################################

prefix_part_names <- c('contrast_variable', 'reference_level', 'treatment_level', 'blocking_variables')
prefix_parts <- unlist(lapply(prefix_part_names, function(x) gsub("[^[:alnum:]]", "_", opt[[x]])))
output_prefix <- paste(prefix_parts[prefix_parts != ''], collapse = '-')

# Load RDS file
dds <- readRDS(opt$deseq_rds)

# Dispersion plot
png(
    file = paste(output_prefix, 'deseq2.dispersion.png', sep = '.'),
    width = opt$plot_width,
    height = opt$plot_height,
    res = opt$plot_res,
    pointsize = opt$plot_point_size
)
plotDispEsts(dds)
dev.off()

# MA plot
png(
    file = paste(output_prefix, 'deseq2.ma.png', sep = '.'),
    width = opt$plot_width,
    height = opt$plot_height,
    res = opt$plot_res,
    pointsize = opt$plot_point_size
)
plotMA(dds)
dev.off()

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
