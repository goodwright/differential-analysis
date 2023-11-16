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

    pca_title = NULL,
    pca_num_genes = 1000,
    loading_top_genes = 10,
    gene_pca_num_genes = 1000,
    gene_pca_alpha = 0.2,
    gene_pca_show_names = FALSE,

    plot_width = 1800,
    plot_height = 1200,
    plot_res = 200,
    plot_point_size = 4
)
opt_types <- lapply(opt, class)


# pairwise correlation gene counts
# genespca top gene count - include default switch for not diisplaying gene names
# choices?
# alpha - use default from other script
# varname.size = 5
# loadings top genes

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

# Check if required parameters have been provided
required_opts <- c('deseq_rds')
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
library(pcaExplorer)

################################################
################################################
## Load data and generate plots               ##
################################################
################################################

# Load RDS file
dds <- readRDS(opt$deseq_rds)

# rlog trans
dds_rlog <- rlogTransformation(dds)

# Create intgroup
intgroup.vars <- c(opt$contrast_variable)
if (!is.null(opt$blocking_variables)) {
    blocking.vars = unlist(strsplit(opt$blocking_variables, split = ';'))
    intgroup.vars <- c(intgroup, blocking.vars)
}

# PCA Plot
pca_title <- ifelse(is.null(opt$pca_title), 'PCA Plot', opt$pca_title)
pdf(
    file = 'pcaexp.pca.pdf'
)
pcaplot(dds_rlog,intgroup = intgroup.vars,ntop = opt$pca_num_genes, pcX = 1, pcY = 2, title = pca_title, ellipse = TRUE)
dev.off()

# Scree plot
pcaobj_dds <- prcomp(t(assay(dds_rlog)))
pdf(
    file = 'pcaexp.scree.pdf'
)
pcascree(pcaobj_dds,type="pev", title=ifelse(is.null(opt$pca_title), 'Scree Plot', opt$pca_title))
dev.off()

# # Corr plot
# dd_corr <- correlatePCs(pcaobj_dds,colData(dds))
# png(
#     file = paste(output_prefix, 'pcaexp.corr.png', sep = '.'),
#     width = opt$plot_width,
#     height = opt$plot_height,
#     res = opt$plot_res,
#     pointsize = opt$plot_point_size
# )
# plotPCcorrs(dd_corr)
# dev.off()

# # extract the table of the genes with high loadings
# hi_loadings(pcaobj_airway,topN = 10,exprTable=counts(dds_airway))

# Plot Loadings
pdf(
    file = 'pcaexp.loadings.pdf'
)
hi_loadings(pcaobj_dds,topN = opt$loading_top_genes)
dev.off()

# Genes PCA
pdf(
    file = 'pcaexp.genes_pca.pdf'
)
genespca(dds_rlog,
        ntop=opt$gene_pca_num_genes,
        alpha = opt$gene_pca_alpha,
        useRownamesAsLabels=opt$gene_pca_show_names,
        varname.size = 3
    )
dev.off()

# Corr plot
pdf(
    file = 'pcaexp.corr.pdf'
)
pair_corr(counts(dds)[1:100,])
dev.off()

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
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
