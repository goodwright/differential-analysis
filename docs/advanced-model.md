# Explicit models, interactions and fixed gene modules

`analysis_mode=design_contrasts` fits gene-level effects with **DESeq2**. Explicitly adding `module_test=camera` and a fixed GMT gene-set file runs a **separate limma-voom CAMERA process** with the same validated design and contrasts. DESeq2 supports paired interactions and numeric contrasts; the previous limitation was the pipeline interface, not DESeq2.

DESeq2 outputs are published under `design_deseq2/`; optional module outputs under `camera_modules/`. CAMERA does not consume DESeq2 p-values or claim its module effects are DESeq2 estimates. Legacy pairwise DESeq2 remains the default, and its plots/GSEA are unchanged. Legacy plots/GSEA are not attached to the new explicit-design mode.

## The paired polysome question

For each biological/fractionation unit, include its input and polysome library. A pair identifier controls shared baseline abundance. Use numeric indicators `poly`, `ko1` and `ko9` (0 or 1):

```
~ 0 + pair + poly + poly:ko1 + poly:ko9
```

Set `numeric_covariates=poly,ko1,ko9`. Other model columns are factors, with sorted levels saved in the audit. Adding genotype main effects to a model whose pair IDs already encode genotype is confounded; the pipeline rejects it. At least two residual degrees of freedom are required. Technical repeats are not automatically independent biological replicates, and a valid design matrix cannot establish biological independence.

Supply explicit coefficient weights as CSV or TSV:

```csv
contrast,coefficient,weight
average,poly:ko1,0.5
average,poly:ko9,0.5
clone1,poly:ko1,1
clone9,poly:ko9,1
```

The average contrast measures the equal-weight mean of `(polysome − input)KO − (polysome − input)WT` across clones. Coefficient names must exactly match R's model matrix. Unsupported expressions, missing values, unknown coefficients, duplicated coefficient rows, all-zero contrasts and rank-deficient models fail rather than silently changing the test. Design formulas permit only column names and `~ + - * : ()`; arbitrary R functions are disallowed. The model and contrast matrices are saved.

## Inputs and modules

Provide a single merged CSV/TSV count matrix with `gene_id`, optional `gene_name`, and exactly the samples in the samplesheet. The first samplesheet columns remain `sample_id,condition` for compatibility with the existing checker; advanced-mode factors follow these. Use raw integer counts or explicitly prepared tximport `lengthScaledTPM` counts with `count_matrix_type=length_scaled_counts`. TPM and log expression are not suitable counts. For DESeq2 only, explicitly declared `length_scaled_counts` estimates are rounded to integers using R `round`, consistent with the tximport-to-DESeq2 count conversion. The audit records the number of fractional entries and maximum rounding change. CAMERA retains the unrounded estimates. Uncorrected Salmon estimated counts need a supported tximport preparation upstream; do not label arbitrary fractional data as lengthScaledTPM. No transcript-length offset is added to counts already prepared with lengthScaledTPM.

Optional GMT `gene_sets` supplies fixed gene modules, not automatically chosen pathways. `gene_set_id_column` must correspond to unique, nonempty IDs in the matrix. Translate the fixed gene list to the reference annotation before outcome inspection; the pipeline does not silently resolve aliases or ambiguous mappings. An optional `gene_universe` text file supplies the prespecified module-plus-background universe **for CAMERA only**. It never restricts the DESeq2 gene-level fit, normalization or multiple-test family. Matching controls by abundance/length/GC is an upstream preparation step, not an automatic feature of CAMERA. Every observed module member must remain in that supplied universe; expression eligibility can still remove unexpressed members.

For an input-controlled protocol, `filter_column=fraction`, `filter_level=input`, `min_cpm=1`, `min_samples=2` controls expression eligibility without using knockout effect direction. This supports the primary paired model but does not by itself enforce a biological protocol's component coverage, clone-agreement or minimum-effect rules; assess those from the saved tables. TMM normalization and voom precision weights use all expression-eligible genes before restricting the competitive universe. This cannot recover absolute global translational changes without spike-ins.

## DESeq2 gene-level tests

The DESeq2 process fits the validated full-rank model matrix once and applies each numeric coefficient contrast through `results()`. It uses Wald tests, `betaPrior=FALSE`, unshrunk log2 fold changes, and normal-approximation 95% confidence intervals. `design_deseq_fit_type` selects `parametric` (default), `local` or `mean`; `design_deseq_sf_type` selects `ratio` (default) or `poscounts`. If DESeq2 automatically changes its dispersion fit, both requested and actual fits are recorded. Legacy `dsq_*` parameters apply to the pairwise path; these explicit-design options are separate.

Normalization and dispersion fitting use all expression-eligible genes. Independent filtering is disabled after the explicit CPM filter, and BH adjustment is applied separately per contrast across nonmissing p-values. Default DESeq2 Cook's-distance filtering is retained: excluded p-values remain NA, with an explicit status and saved Cook's distances. Automatic count replacement is disabled (`minReplicatesForReplace=Inf`). A nonconverged coefficient fit stops the process rather than reporting an apparently complete result. A saved coefficient-name mapping connects the input contrast names to DESeq2's internal names.

## Separate fixed-module test

Set `module_test=camera` explicitly; the default `none` runs DESeq2 alone. A gene-set file without this setting, or this setting without a gene-set file, is rejected to avoid silently choosing a statistical test. CAMERA fits the same design and contrasts to TMM-normalized voom expression in its own process. Its null hypothesis is competitive: module genes do not shift more than genes in the declared background. This differs from the individual-gene DESeq2 null. Gene-set membership and background must be fixed before looking at validation outcomes.

CAMERA estimates intergene correlation with `inter.gene.cor=NA`, using finite residual degrees of freedom for the estimated-correlation test. Negative estimated correlations are truncated for the variance-inflation calculation; the original estimates remain in the output. With observation-specific voom weights, CAMERA estimates can vary slightly under different contrast parameterizations, so module p-values for a contrast and its negative need not be numerically identical. Gene-level two-sided p-values are invariant to contrast sign. Do not plug an estimated correlation into CAMERA's fixed-correlation mode: that changes its degrees of freedom and overstates precision in small experiments. CAMERA reports direction and p-value. Module p-values are adjusted across every tested module × contrast in that execution using the requested `BH`, `holm` or `bonferroni` method. DESeq2 gene-level correction is a separate family and is never pooled with module p-values. For one prespecified primary module/contrast, submit that as its own test family, and place secondary contrasts in a separate declared family. Do not choose the family after seeing outcomes.

CAMERA effect summaries use its own limma-voom estimates, labeled in the output, and include module and background mean log2 changes, their difference, and the fraction of module genes with negative changes. These summaries have no module-effect confidence interval; gene-level confidence intervals and the correlation-aware CAMERA p-value are separate outputs. A low p-value alone is not the biological decision rule. Missing/filtered members and untested sets are explicitly reported.

## Runtime and Flow registration

The new mode uses DESeq2 1.44.0, limma 3.60.6 and edgeR 4.2.2 with R 4.4.1 (Bioconductor 3.19). The Dockerfile pins its Bioconductor base by digest and asserts these versions. Every built image records all installed packages in `/opt/advanced-model/installed-packages.tsv` and R session information in `/opt/advanced-model/sessionInfo.txt`. Dependency repositories can change between builds; the published image digest, rather than a rebuild of the Dockerfile, identifies the exact runtime.

The Advanced runtime container workflow builds for Linux amd64, runs the statistical regression suite inside the image, then runs both the DESeq2-plus-CAMERA and DESeq2-only Nextflow paths using Docker. It compares each pipeline result with the direct R result produced in that same image. Only after these tests pass can a same-repository PR publish a commit-tagged candidate to GHCR; fork PRs build and test without publishing. Publication does not deploy the pipeline or change a mutable production tag. Build logs and uploaded runtime evidence preserve package versions, the image digest and test results. New GHCR packages can be private: verify anonymous pull access before using a digest in Flow.

The default image is `ghcr.io/goodwright/differential-analysis/advanced-model@sha256:e272a5038cc3a80acf1cefc826658fbb6b48efc60ee1f7729affeda583d8db14`, built and tested by [container run 37043635091](https://github.com/goodwright/differential-analysis/actions/runs/37043635091) from source `bdf59b435aaa973351d00db5e39f3ee47b45d4f7`. Both Nextflow Docker paths matched their direct R reference results. **Deployment prerequisite:** this new GHCR package was initially private, and an anonymous pull check returned HTTP 401. A package administrator must make it public, or Flow's execution environment must have separately configured pull access, before deployment. A successful build and publication do not establish Flow access.

`flow/schema/main.json` is the single canonical Flow registration schema. It uses an input-section array, `params`, uploaded-file inputs of type `data`, custom samplesheet columns and output declarations. The duplicate registration candidate has been removed to prevent the wrong format being deployed again. A schema contract test guards this format and the registered runtime default. Register a new Flow version linked to the reviewed source commit and tested immutable container; avoid overwriting existing releases. Registration and biological execution are separate from this code change.

Local tests can use an installed R environment with the exact package versions above. Conda execution is rejected for the advanced mode because these exact Bioconductor patch versions are unavailable in the channels; existing pairwise Conda support is unchanged.

## Tests

`tests/advanced_model/test_advanced.R` creates synthetic paired libraries, simulates a module-specific polysome effect and a separate total-RNA effect, and checks contrast direction, nuisance separation, sample order, correction families, missing-module handling and invalid designs. It never reads the held-out GEO validation data.
