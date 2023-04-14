# differential-analysis

## Pipeline Summary

Differential analysis is a statistical method used to identify differences in gene expression or other molecular features between two or more groups of samples. This approach is critical in identifying molecular biomarkers associated with disease or other biological conditions, as well as understanding the underlying mechanisms of complex biological processes.

This pipeline is designed to process simple count matrices, which are the output of many count-based workflows such as RNA-seq. Count matrices for other applications can also be easily created as a simple sample versus feature count table and used as input. The count matrix input is paired with a sample sheet that describes metadata associated with each sample. The sample sheet can include any number of attributes that group the sample experimentally.

A set of comparisons based on the supplied statistical design, using reference and contrast factors, as well as more complex comparisons involving blocking factors is then performed. The primary output is a results table for each comparison, which shows the statistical significance of each feature. We also provide a rich set of reporting outputs to help users assess the quality of their data and visualise their results.

The pipeline executes the following main stages:

- Input sample sheet checking and experimental validation
- Run differential analysis using DESeq2
- Generate DESeq2 standard plots
- Run dimenionality reduction on the DESeq2 output using pcaExplorer
- Run geneset enrichment analysis on each comparison result set

## Quick-start

To test the pipeline run it with the `test` profile and the container engine you wish to use eg. `docker`. For example:

```bash
nextflow run main.nf -profile test,docker
```

## Tool Usage List

- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) - Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.
- [pcaExplorer](https://www.bioconductor.org/packages/release/bioc/html/pcaExplorer.html) - Functionality for visualization of RNA-seq datasets based on Principal Components Analysis.
- [genekitr](https://www.genekitr.fun/) - A collection of R tools designed to aid in performing and visualising gene set enrichment analysis

## Input

### Count Matrix

Supplied using the `--counts` parameter. This is a numeric square matrix file, comma or tab-separated, with a column for every sample against a count of features. For example:

| gene_id            | gene_name | SRX8042381 | SRX8042382 | SRX8042383 |
| ------------------ | --------- | ---------- | ---------- | ---------- |
| ENSMUSG00000023978 | Prph2     | 10         | 2          | 4          |
| ENSMUSG00000030324 | Rho       | 181        | 14         | 12         |
| ENSMUSG00000031450 | Grk1      | 20         | 3          | 1          |
| ENSMUSG00000056055 | Sag       | 52         | 3          | 10         |
| ENSMUSG00000025900 | Rp1       | 17         | 3          | 2          |
| ENSMUSG00000026834 | Acvr1c    | 10         | 4          | 366        |
| ENSMUSG00000025386 | Pde6g     | 7          | 8          | 1          |

At a minmimum an id column must be provided which is always defined as the first column. Additional info columns can be included for context which can be utilised by downstream processes such as gene name. Meta data is seperated from sample columns by the sample ids provided in the `samplesheet`

### Samplesheet

Supplied using the `--samplesheet` parameter. At a minimum the samplesheet must contain a set of sample ids which match what is in the supplied counts file and at least one column to use as a comparison variable. Any number of additional status columns can be included for more complex experimental designs such as in the example below where an experimental condition and batch number are provided.

| sample_id           | condition | batch |
| ------------------- | --------- | ----- |
| WT_REP1             | A         | 1     |
| WT_REP2             | A         | 2     |
| RAP1_UNINDUCED_REP1 | B         | 1     |
| RAP1_UNINDUCED_REP2 | B         | 2     |

### Experimental Design

The experimental design that defines which comparisons and blocking factors are used can be tuned using three different parameters:

- `contrast_column`: This is the column in the samplesheet from which comparisons will be generated. By default this is set to `condition`.
- `comparisons`: This defines which comparisons will be run from the contrast column. By default this is set to `all`, in this configuration an all against all comparison set of all variables in the contrast column will be performed. In designs that have a large number of variables or where some comparisons would be invalid, users can defined their own comparison set using the following notation `C_A:B_C`. Comparisons must be delimited by `:` with the reference factor on the left seperated by an underscore with the contrast value on the right.
- `blocking_factors`: This variable allows for multi-factor experimental designs. By default this is blank. Blocking factors can be defined by a `;` seperated list of factors that exist in the input samplesheet. See the DESeq2 docs on multi-factor design [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs) for more information.

## Output

Pipeline output will default to the path `./results` unless modified with the `--outdir` parameter.

The `deseq2` folder contains output data pertaining to the differential analysis comparisons. The folder will always contain the size factors and normalised counts in its root folder as well as a `plots` folder that contains experiment-level diagnostic plots such as the dispersion and ma plots. 

For each comparison an additional sub-folder will be included that contains the primary output results table and any other comparison specific plots such as Volcano plots.

Additional folders are outputted depending on which additional processing modules are active and are clearly named. For example the geneset enrichment analysis module is names `gsea`.

## Pipeline Parameters

- `--study_name`: All output from the pipeline will be placed in the directory `$outdir/$study_name`. This is so several runs with different parameters can be performed without needing to worry about mixing up results. The default path is `results/study`.

- `--count_sep`: Defines the count separator used in the input count files. Defaults to tab `\t`

### Differential Analysis

This pipeline exposes most key DESeq2 parameters via nextflow parameters with the `dsq_` prefix. More information on these parameters can be found in the table below, or in the DESeq2 documentation.

| Parameter                      | Default    | Description                                                                           |
| ------------------------------ | ---------- | ------------------------------------------------------------------------------------- |
| dsq_test                       | Wald       | Options: Wald, LRT                                                                    |
| dsq_fit_type                   | parametric | Options: parametric, local, mean, glmGamPoi                                           |
| dsq_min_replicates_for_replace | 7          | The minimum number of replicates required in order to use replaceOutliers on a sample |
| dsq_sf_type                    | ratio      | Options: ratio, poscounts, iterate                                                    |
| dsq_lfc_threshold              | 0          | A non-negative value which specifies a log2 fold change threshold                     |
| dsq_alt_hypothesis             | greaterAbs | Options: greaterAbs, lessAbs, greater, less                                           |
| dsq_alpha                      | 0.1        | The significance cutoff used for optimizing the independent filtering                 |
| dsq_minmu                      | 0.5        | Lower bound on the estimated count (used when calculating contrasts)                  |
| dsq_lfcshrink_type             | ashr       | Options: apeglm, ashr, normal                                                         |
| dsq_p_thresh                   | 1          | p-value cut off for filtering the deseq2 results table for downstream modules         |

### Gene set Enrichment Analysis

The most important parameter for this module is correct setting of `--organism`, this informs genekitr which organism to pull data for and must be set by finding the correct species name [here](https://genekitr.online/docs/species.html). If no organism is suitable or users wish to switch this module off then the `--skip_gsea` parameter can be used. genekitr parameters are exposed via nextflow parameters with the `gsea_` prefix. More information on these parameters can be found in the table below, or in the genekitr documentation.

| Parameter             | Default    | Description                                                                       |
| --------------------- | ---------- | --------------------------------------------------------------------------------- | ------------------------ | ------------------------ |
| gsea_p_cutoff         | 0.05       | A numeric of cutoff for both pvalue and adjusted pvalue                           |
| gsea_q_cutoff         | 0.05       | A numeric of cutoff for both qvalue                                               |
| gsea_ontology         | mf         | Biological Processes (BP)                                                         | Molecular Functions (MF) | Cellular Components (CC) |
| gsea_min_gset_size    | 10         | Minimal size of each gene set for analysis                                        |
| gsea_max_gset_size    | 500        | Max size of each gene set for analysis                                            |
| gsea_p_adjust_method  | BH         | Choose from “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none” |
| gsea_pathway_count    | 10         | How many pathways to show on the plots                                            |
| gsea_stats_metric     | p.adjust   | Stats metric for the plots - p.adjust, pvalue, qvalue                             |
| gsea_term_metric      | FoldEnrich | Term metric for the ora plots - FoldEnrich, GeneRatio, Count, RichFactor          |
| gsea_scale_ratio      | 0.25       | Plot scale ratio                                                                  |
| gsea_main_text_size   | 5          | Plot text size                                                                    |
| gsea_legend_text_size | 8          | Plot legend size                                                                  |

## Authors and contact

This DSL2 Nextflow pipeline is maintained by Goodwright.
To raise any issues or comments with the pipeline you can (in order of preference):

- Raise an issue in this repository
- Write to us in our [Slack](https://join.slack.com/t/imapsgroup/shared_invite/zt-r24y3591-Xbhnym2t38u_urU~I0K0lQ)
- Email charlotte.capitanchik@goodwright.com
