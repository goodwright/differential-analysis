#!/usr/bin/env bash
# Run both public workflow paths inside the candidate image before publishing it.
set -euo pipefail
image_ref="$1"
fixture_dir="$2"
output_dir="$3"
for mode in camera none; do
    module_args=()
    if [[ "$mode" == camera ]]; then
        module_args=(--module_test camera --gene_sets "$fixture_dir/modules.gmt" --gene_universe "$fixture_dir/universe.txt")
    fi
    nextflow run main.nf -profile docker \
        --analysis_mode design_contrasts --custom_model_container "$image_ref" \
        --samplesheet "$fixture_dir/samples.csv" --counts "$fixture_dir/counts.tsv" \
        --design_formula '~ 0 + pair + poly + poly:ko1 + poly:ko9' \
        --numeric_covariates poly,ko1,ko9 --contrast_table "$fixture_dir/contrasts.csv" \
        "${module_args[@]}" --filter_column fraction --filter_level input \
        --study_name "container_$mode" --outdir "$output_dir" \
        --publish_dir_mode copy --max_memory 4.GB --max_cpus 2
    diff "$fixture_dir/deseq_result/average.deseq2.results.tsv" \
        "$output_dir/container_$mode/design_deseq2/average.deseq2.results.tsv"
    if [[ "$mode" == camera ]]; then
        diff "$fixture_dir/result/camera.results.tsv" \
            "$output_dir/container_$mode/camera_modules/camera.results.tsv"
    else
        test ! -d "$output_dir/container_$mode/camera_modules"
    fi
done
