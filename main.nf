#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    goodwright/differential_analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/goodwright/differential_analysis
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { summary_log } from './modules/goodwright/util/logging/main'
// include { multiqc_summary } from './modules/goodwright/util/logging/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info summary_log(workflow, params, params.debug, params.monochrome_logs)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    samplesheet: params.samplesheet,
    counts     : params.counts
]
for (param in check_param_list) { 
    if (!param.value) { 
        exit 1, "Required parameter not specified: ${param.key}"
    } 
    else {
        file(param.value, checkIfExists: true)
    }
}

// // Check non-manditory input parameters to see if the files exist if they have been specified
// check_param_list = [
//     params.target_genome_index,
//     params.smrna_genome_index
// ]
// for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Collect comparisons if any specified
comparisons = params.comparisons ? params.comparisons.split(':').collect{ it.trim() } : null

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

// include { MULTIQC } from './modules/local/multiqc'

//
// SUBWORKFLOWS
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT GOODWRIGHT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

include { SAMPLE_DIFF_SAMPLESHEET_CHECK } from './modules/goodwright/sample/diff_samplesheet_check/main'
include { R_DESEQ2                      } from './modules/goodwright/r/deseq2/main'
include { R_DESEQ2_PLOTS                } from './modules/goodwright/r/deseq2_plots/main'
include { R_PCAEXPLORER                 } from './modules/goodwright/r/pcaexplorer/main'
// include { DUMP_SOFTWARE_VERSIONS } from './modules/goodwright/dump_software_versions/main'

//
// SUBWORKFLOWS
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

//
// SUBWORKFLOWS
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIFF_ANALYSIS {
    // Init
    ch_versions = Channel.empty()
    // ch_target_genome_index = []
    // ch_smrna_genome_index  = []

    // Prepare manditory params into channels 
    ch_samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_counts      = file(params.counts, checkIfExists: true)

    // // Prepare non-manditory params into channels
    // if(params.target_genome_index) { ch_target_genome_index = file(params.target_genome_index, checkIfExists: true) }
    // if(params.smrna_genome_index)  { ch_smrna_genome_index = file(params.smrna_genome_index, checkIfExists: true) }

    ch_meta = Channel.empty()
    if(params.run_input_check) {
        /*
        * MODULE: Check the samplesheet for errors
        */
        SAMPLE_DIFF_SAMPLESHEET_CHECK (
            ch_samplesheet,
            ch_counts,
            params.count_sep
        )
        ch_versions = ch_versions.mix(SAMPLE_DIFF_SAMPLESHEET_CHECK.out.versions)

        /*
        * MODULE: Parse samplesheet into meta and fastq files
        */
        ch_meta = SAMPLE_DIFF_SAMPLESHEET_CHECK.out.csv
           .splitCsv ( header:true, sep:"," )
    }
    //EXAMPLE CHANNEL STRUCT: [sample_id:RAP1_IAA_30M_REP1, condition:B]
    //ch_meta | view
    //SAMPLE_DIFF_SAMPLESHEET_CHECK.out.csv | view

    if(params.run_diff_analysis) {
        /*
        * CHANNEL: Create channel from samplesheet
        */
        ch_design = SAMPLE_DIFF_SAMPLESHEET_CHECK.out.csv
            .collect()
            .map { ["", it]}
        //ch_design | view

        /*
        * CHANNEL: Create channel for all against all analysis
        *         but filter for only the conditions specified
        */
        ch_comparison_set = ch_meta
            .map { it[params.contrast_column] }
            .unique()
        //ch_comparison_set | view

        ch_comparisons = ch_comparison_set
            .combine(ch_comparison_set)
            .filter { it[0] != it[1] }
        //ch_comparisons | view

        if( comparisons ) {
            ch_comparisons = ch_comparisons
                .filter { ( it[0] + "_" + it[1] )  in comparisons }
        }
        //ch_comparisons | view

        /*
        * MODULE: Run deseq2
        */
        R_DESEQ2 (
            ch_design.collect(),
            ch_counts,
            params.contrast_column,
            ch_comparisons.map { it[0] },
            ch_comparisons.map { it[1] },
            params.blocking_factors
        )
        ch_versions = ch_versions.mix(R_DESEQ2.out.versions)

        /*
        * MODULE: Run deseq2 plots
        */
        R_DESEQ2_PLOTS (
            R_DESEQ2.out.rdata,
            params.contrast_column,
            ch_comparisons.map { it[0] },
            ch_comparisons.map { it[1] },
            params.blocking_factors
        )
        ch_versions = ch_versions.mix(R_DESEQ2_PLOTS.out.versions)
    }



    // if(params.run_reporting) {
    //     /*
    //     * MODULE: Collect software versions
    //     */
    //     DUMP_SOFTWARE_VERSIONS (
    //         ch_versions.unique().collectFile()
    //     )

    //     /*
    //     * MODULE: Run multiqc
    //     */
    //     workflow_summary    = multiqc_summary(workflow, params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    //     MULTIQC (
    //         ch_multiqc_config,
    //         DUMP_SOFTWARE_VERSIONS.out.mqc_yml.collect(),
    //         DUMP_SOFTWARE_VERSIONS.out.mqc_unique_yml.collect(),
    //         ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yml"),
    //         FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
    //         FASTQC_TRIMGALORE.out.fastqc_trim_zip.collect{it[1]}.ifEmpty([]),
    //         FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
    //         ch_bt_log.collect{it[1]}.ifEmpty([]),
    //         ch_star_log.collect{it[1]}.ifEmpty([]),
    //         CLIPSEQ_CLIPQC.out.tsv.collect().ifEmpty([])
    //     )
    // }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    DIFF_ANALYSIS ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EVENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/