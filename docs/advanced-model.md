# Explicit models, interactions and fixed gene modules

`analysis_mode=design_contrasts` adds a limma-voom model and optional CAMERA competitive gene-set testing. Pairwise DESeq2 remains the default mode. This mode is a separate statistical engine, not DESeq2 interaction coefficients relabeled as limma results. Legacy DESeq2 plotting/GSEA modules are not run for the advanced mode; its result and audit tables are published under `custom_model/`.

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

Provide a single merged CSV/TSV count matrix with `gene_id`, optional `gene_name`, and exactly the samples in the samplesheet. The first samplesheet columns remain `sample_id,condition` for compatibility with the existing checker; advanced-mode factors follow these. Use raw integer counts or explicitly prepared tximport `lengthScaledTPM` counts with `count_matrix_type=length_scaled_counts`. TPM and log expression are not suitable counts. Fractional counts are never silently rounded in advanced mode.

Optional GMT `gene_sets` supplies fixed gene modules, not automatically chosen pathways. `gene_set_id_column` must correspond to unique, nonempty IDs in the matrix. Translate the fixed gene list to the reference annotation before outcome inspection; the pipeline does not silently resolve aliases or ambiguous mappings. An optional `gene_universe` text file supplies the prespecified module-plus-background universe. Matching controls by abundance/length/GC is an upstream preparation step, not an automatic feature of CAMERA. Every observed module member must remain in that supplied universe; expression eligibility can still remove unexpressed members.

For an input-controlled protocol, `filter_column=fraction`, `filter_level=input`, `min_cpm=1`, `min_samples=2` controls expression eligibility without using knockout effect direction. This supports the primary paired model but does not by itself enforce a biological protocol's component coverage, clone-agreement or minimum-effect rules; assess those from the saved tables. TMM normalization and voom precision weights use all expression-eligible genes before restricting the competitive universe. This cannot recover absolute global translational changes without spike-ins.

CAMERA estimates intergene correlation with `inter.gene.cor=NA`, using finite residual degrees of freedom for the estimated-correlation test. Negative estimated correlations are truncated for the variance-inflation calculation; the original estimates remain in the output. With observation-specific voom weights, CAMERA estimates can vary slightly under different contrast parameterizations, so module p-values for a contrast and its negative need not be numerically identical. Gene-level two-sided p-values are invariant to contrast sign. Do not plug an estimated correlation into CAMERA's fixed-correlation mode: that changes its degrees of freedom and overstates precision in small experiments. CAMERA reports direction and p-value. Module p-values are adjusted across every tested module × contrast in that execution using the requested `BH`, `holm` or `bonferroni` method. Gene-level BH correction is performed separately per contrast. For one prespecified primary module/contrast, submit that as its own test family, and place secondary contrasts in a separate declared family. Do not choose the family after seeing outcomes.

Effect summaries include module and background mean log2 changes, their difference, and the fraction of module genes with negative changes. These summaries have no module-effect confidence interval; gene-level confidence intervals and the correlation-aware CAMERA p-value are separate outputs. A low p-value alone is not the biological decision rule. Missing/filtered members and untested sets are explicitly reported.

## Runtime and Flow registration

The new mode targets limma 3.60.6 / edgeR 4.2.2 with R 4.4.1. An environment specification and Docker build recipe live in `lib/advanced_model/`. The image must be built, tested and pinned by digest in `custom_model_container` before a production Flow release. The Dockerfile's base tag is not an immutable release artifact; capture its resolved digest and all installed dependency versions when building. Local tests can use a pinned installed R environment.

`flow/schema/main.json` retains the repository's legacy schema format and includes the added fields. `flow/schema/advanced-flow-schema.json` extends the live Flow schema retrieved on 2026-10-02, including custom sample columns and output declarations. This is a registration candidate, not evidence that the live pipeline has been updated. Register a new version linked to the reviewed source commit and immutable container; never overwrite the existing 1.0 release. No registration, deployment or biological execution is part of this code change.

## Tests

`tests/advanced_model/test_advanced.R` creates synthetic paired libraries, simulates a module-specific polysome effect and a separate total-RNA effect, and checks contrast direction, nuisance separation, sample order, correction families, missing-module handling and invalid designs. It never reads the held-out GEO validation data.
