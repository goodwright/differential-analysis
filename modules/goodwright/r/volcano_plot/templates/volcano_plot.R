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
        row.names = row.names
    )
}

#####################################################
#####################################################
## PARSE PARAMETERS FROM NEXTFLOW AND COMMAND LINE ##
#####################################################
#####################################################

# Set defaults and classes
opt <- list(
    deseq_results = "!{deseq_results}",
    contrast_variable = "!{contrast_variable}",
    reference_level = "!{reference_level}",
    treatment_level = "!{treatment_level}",
    blocking_variables = "!{blocking_variables}",

    plot_width = 1800,
    plot_height = 1200,
    plot_res = 300,

    fold_change = !{fold_change},
    p_value = !{p_value}
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
required_opts <- c('deseq_results', 'contrast_variable', 'reference_level', 'treatment_level')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('deseq_results')){
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
library(ggplot2)
library(ggrepel)
library(dplyr)

################################################
################################################
## Load data and generate plots               ##
################################################
################################################

prefix_part_names <- c('contrast_variable', 'reference_level', 'treatment_level', 'blocking_variables')
prefix_parts <- unlist(lapply(prefix_part_names, function(x) gsub("[^[:alnum:]]", "_", opt[[x]])))
output_prefix <- paste(prefix_parts[prefix_parts != ''], collapse = '-')
contrast.name <- paste(opt$treatment_level, opt$reference_level, sep = "_vs_")

# Load results table
de <- read_delim_flexible(
    file = opt$deseq_results,
    header = TRUE
)

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# get numbers of diff genes for labelling
n_unchanged <- de %>% filter(log2FoldChange < opt$fold_change & log2FoldChange > -(opt$fold_change)) %>% nrow()
n_up <- de %>% filter(log2FoldChange >= opt$fold_change & padj < opt$p_value) %>% nrow()
n_down <- de %>% filter(log2FoldChange <= -(opt$fold_change) & padj < opt$p_value) %>% nrow()

de$diffexpressed <- paste0("Unchanged (",n_unchanged,")")
de$diffexpressed[de$log2FoldChange >= opt$fold_change & de$padj < opt$p_value] <- paste0("Up (",n_up,")")
de$diffexpressed[de$log2FoldChange <= -(opt$fold_change) & de$padj < opt$p_value] <- paste0("Down (",n_down,")")

# set colours vector
if (n_up == 0 & n_down == 0){
    cvec = c("#84A1AB")
} else if (n_up == 0){
    cvec = c("#B02302", "#84A1AB")
} else if (n_down == 0){
    cvec = c("#84A1AB", "#61B002")
} else {
    cvec = c("#B02302", "#84A1AB", "#61B002")
}

# label genes that are differentially expressed
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_id[de$diffexpressed != "NO"]

# Volcano plot
ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), label=delabel)) +
        geom_vline(xintercept=c(-(opt$fold_change), opt$fold_change), col="light grey", linetype="dashed") +
        geom_hline(yintercept=-log10(opt$p_value), col="light grey", linetype="dashed") +
        geom_point(aes(color=diffexpressed), alpha=0.5) +
        geom_label_repel(size=3) +
        scale_color_manual(values=cvec) +
        theme_bw()

ggsave(
    file = paste(output_prefix, 'deseq2.volcano.pdf', sep = '.')
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
