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
    cores = 1,                            # Number of cores to use
    deseq2_results = "!{deseq2_results}", # The input deseq2 results table
    prefix = "!{prefix}",                 # Output prefix

    # Results params
    dsq_p_thresh = 1, # Filter threshold for the deseq2 results table

    # General GSEA params
    gsea_p_cutoff = 0.05, # numeric of cutoff for both pvalue and adjusted pvalue, default should be 0.05
    gsea_q_cutoff = 0.05, # numeric of cutoff for qvalue, default should be 0.15

    # genekitr params
    ontology = "mf",            # Biological Processes (BP) | Molecular Functions (MF) | Cellular Components (CC)
    organism = "!{organism}",   # https://genekitr.online/docs/species.html
    min_gset_size = 10,         # Minimal size of each gene set for analysis. Default should be 10
    max_gset_size = 500,        # Maximal size. Default should be 500.
    p_adjust_method = "BH",     # Choose from “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”
    pathway_count = 10,         # How many pathways to show on the plots
    stats_metric = "p.adjust",  # Stats metric for the plots - c("p.adjust", "pvalue", "qvalue")
    term_metric = "FoldEnrich", # Term metric for the ora plots - c("FoldEnrich", "GeneRatio", "Count", "RichFactor")
    scale_ratio = 0.25,         # Plot scale ratio
    main_text_size = 5,         # Plot text size
    legend_text_size = 8,       # Plot legend size

    # General Plotting params
    plot_width = 1800,
    plot_height = 1200,
    plot_res = 300
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

# Set nulls for required params that havent been resolved
if (startsWith(opt$organism, "!")) { opt$organism <- NULL }
if (startsWith(opt$prefix, "!")) { opt$prefix <- NULL }

# Check if required parameters have been provided
required_opts <- c('deseq2_results', 'organism', 'prefix')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]

if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid
for (file_input in c('deseq2_results')){
    if (is.null(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }

    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

# Convert params
opt$dsq_p_thresh <- as.numeric(opt$dsq_p_thresh)
opt$gsea_p_cutoff <- as.numeric(opt$gsea_p_cutoff)
opt$gsea_q_cutoff <- as.numeric(opt$gsea_q_cutoff)
opt$min_gset_size <- as.numeric(opt$min_gset_size)
opt$max_gset_size <- as.numeric(opt$max_gset_size)
opt$pathway_count <- as.numeric(opt$pathway_count)
opt$scale_ratio <- as.numeric(opt$scale_ratio)
opt$main_text_size <- as.numeric(opt$main_text_size)
opt$legend_text_size <- as.numeric(opt$legend_text_size)

print(opt)

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(geneset)
library(genekitr)
library(patchwork)
library(igraph)
library(ggraph)

################################################
################################################
## Load and prepare data                      ##
################################################
################################################

# Load results table
results <- read_delim_flexible(
    file = opt$deseq2_results,
    header = TRUE
)

# Filter table based on p-value
results <- subset(results, padj < opt$dsq_p_thresh)

# summary(results)

# Sort by descending logfold change
results <- results[order(results$log2FoldChange, decreasing = TRUE),]

# Extract named vector of genes and fold change
gene_list <- results$log2FoldChange
names(gene_list) <- results$gene_id

# Get up genes
results_up <- subset(results, log2FoldChange >= 0)

#Get down genes
results_down <- subset(results, log2FoldChange < 0)

################################################
################################################
## Generate genekitr                          ##
################################################
################################################

gs <- geneset::getGO(org=opt$organism, ont=opt$ontology)

gse <- genGSEA(
    genelist=gene_list,
    geneset=gs,
    padj_method=opt$p_adjust_method,
    min_gset_size=opt$min_gset_size,
    max_gset_size=opt$max_gset_size,
    p_cutoff=opt$gsea_p_cutoff,
    q_cutoff=opt$gsea_q_cutoff
)

ora <- genORA(
    results$gene_id,
    geneset=gs,
    padj_method=opt$p_adjust_method,
    min_gset_size=opt$min_gset_size,
    max_gset_size=opt$max_gset_size,
    p_cutoff=opt$gsea_p_cutoff,
    q_cutoff=opt$gsea_q_cutoff
)

ora_up <- genORA(
    results_up$gene_id,
    geneset=gs,
    padj_method=opt$p_adjust_method,
    min_gset_size=opt$min_gset_size,
    max_gset_size=opt$max_gset_size,
    p_cutoff=opt$gsea_p_cutoff,
    q_cutoff=opt$gsea_q_cutoff
)

ora_down <- genORA(
    results_down$gene_id,
    geneset=gs,
    padj_method=opt$p_adjust_method,
    min_gset_size=opt$min_gset_size,
    max_gset_size=opt$max_gset_size,
    p_cutoff=opt$gsea_p_cutoff,
    q_cutoff=opt$gsea_q_cutoff
)

ora_filt <- head(ora, opt$pathway_count)
ora_up_filt <- head(ora_up, opt$pathway_count)
ora_down_filt <- head(ora_down, opt$pathway_count)

################################################
################################################
## Output Data                                ##
################################################
################################################

genekitr::expoSheet(data_list = gse,
                    data_name = names(gse),
                    filename = paste(opt$prefix, "genekitr_gsea_result.xlsx", sep = '.'),
                    dir = "./")

genekitr::expoSheet(data_list = ora,
                    data_name = names(ora),
                    filename = paste(opt$prefix, "genekitr_ora_result.xlsx", sep = '.'),
                    dir = "./")

################################################
################################################
## Output GSEA Plots                          ##
################################################
################################################

gse$gsea_df[,1] <- gse$gsea_df$Description
output_count <- length(gse$gsea_df[,1])
print("GSE pathway result count")
print(output_count)

if(output_count == 0) {
    print("No GSE pathways detected")
} else {
    # Set pathway count to available ids/2
    up_down_pathway_count <-opt$pathway_count
    if(floor(output_count / 2) < opt$pathway_count) {
        up_down_pathway_count <- floor(output_count / 2)
        print("Adjusted pathway count")
        print(opt$pathway_count)
    }

    # Volcano plot
    png(
        file = paste(opt$prefix, 'gsea.volcano.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotGSEA(gse, plot_type = "volcano", show_pathway = up_down_pathway_count, stats_metric = opt$stats_metric)
    dev.off()

    # Multi-pathway plot
    png(
        file = paste(opt$prefix, 'gsea.mpathway.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotGSEA(gse, plot_type = "fgsea", show_pathway = up_down_pathway_count, stats_metric = opt$stats_metric)
    dev.off()

    # Ridge plot
    png(
        file = paste(opt$prefix, 'gsea.ridge.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotGSEA(gse, plot_type = "ridge", show_pathway = opt$pathway_count, stats_metric = opt$stats_metric)
    dev.off()

    # Bar plot
    png(
        file = paste(opt$prefix, 'gsea.bar.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotGSEA(gse, plot_type = "bar", show_pathway = opt$pathway_count, stats_metric = opt$stats_metric)
    dev.off()
}

################################################
################################################
## Output ORA Plots                          ##
################################################
################################################

output_count <- length(ora[,1])
print("ORA pathway result count")
print(output_count)

if(output_count == 0) {
    print("No ORA pathways detected")
} else {
    # Bar
    png(
        file = paste(opt$prefix, 'ora.bar.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(ora_filt, plot_type = "bar", stats_metric = opt$stats_metric, term_metric = opt$term_metric)
    dev.off()

    # Bubble
    png(
        file = paste(opt$prefix, 'ora.bubble.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(
        ora_filt,
        plot_type = "bubble",
        stats_metric = opt$stats_metric,
        term_metric = opt$term_metric,
        scale_ratio = opt$scale_ratio,
        main_text_size = opt$main_text_size,
        legend_text_size = opt$legend_text_size
    )
    dev.off()

    # Dot
    png(
        file = paste(opt$prefix, 'ora.dot.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(
        ora_filt,
        plot_type = "dot",
        stats_metric = opt$stats_metric,
        term_metric = opt$term_metric,
        scale_ratio = opt$scale_ratio,
        main_text_size = opt$main_text_size,
        legend_text_size = opt$legend_text_size
    )
    dev.off()

    # Lollipop
    png(
        file = paste(opt$prefix, 'ora.lollipop.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(
        ora_filt,
        plot_type = "lollipop",
        stats_metric = opt$stats_metric,
        term_metric = opt$term_metric,
        scale_ratio = opt$scale_ratio,
        main_text_size = opt$main_text_size,
        legend_text_size = opt$legend_text_size
    )
    dev.off()

    # Geneheat
    png(
        file = paste(opt$prefix, 'ora.geneheat.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(
        ora_filt,
        plot_type = "geneheat",
        stats_metric = opt$stats_metric,
        term_metric = opt$term_metric,
        scale_ratio = opt$scale_ratio,
        main_text_size = opt$main_text_size,
        legend_text_size = opt$legend_text_size
    )
    dev.off()

    # Network
    png(
        file = paste(opt$prefix, 'ora.network.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(
        ora_filt,
        plot_type = "network",
        stats_metric = opt$stats_metric,
        term_metric = opt$term_metric,
        scale_ratio = opt$scale_ratio,
        main_text_size = opt$main_text_size,
        legend_text_size = opt$legend_text_size
    )
    dev.off()

    # Gomap
    png(
        file = paste(opt$prefix, 'ora.gomap.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(
        ora_filt,
        plot_type = "gomap",
        stats_metric = opt$stats_metric,
        term_metric = opt$term_metric
    )
    dev.off()

    # Goheat
    png(
        file = paste(opt$prefix, 'ora.goheat.png', sep = '.'),
        width = opt$plot_width,
        height = opt$plot_height,
        res = opt$plot_res
    )
    plotEnrich(
        ora_filt,
        plot_type = "goheat",
        stats_metric = opt$stats_metric,
        term_metric = opt$term_metric
    )
    dev.off()
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(opt$prefix, "R_sessionInfo.log", sep = '.'))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
genekitr.version <- as.character(packageVersion('genekitr'))

writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', r.version),
        paste('    genekitr:', genekitr.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
