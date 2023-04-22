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
    samplesheet: params.samplesheet
]
for (param in check_param_list) {
    if (!param.value) {
        exit 1, "Required parameter not specified: ${param.key}"
    }
    else {
        file(param.value, checkIfExists: true)
    }
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Collect comparisons if any specified
comparisons = params.comparisons ? params.comparisons.split(':').collect{ it.trim() } : null

// Split count file input into list
count_files = params.counts.split(',').collect{ file(it.trim(), checkIfExists: true) }.flatten()

// Collect blocking variable
ch_blocking_factors = params.blocking_factors ? params.blocking_factors.split(',').collect{ it.trim() } : null

// Parse blocking factors into a channel
if (ch_blocking_factors) {
    ch_blocking_factors = Channel.from(ch_blocking_factors)

    ch_blocking_factors = ch_blocking_factors.map{
        def bsplit = it.split(':')
        [bsplit[0], bsplit.size() == 1 ? '' : bsplit[1]]
    }
}
//ch_blocking_factors | view

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
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
include { R_VOLCANO_PLOT                } from './modules/goodwright/r/volcano_plot/main'
include { R_GSEA                        } from './modules/goodwright/r/gsea/main'
include { DUMP_SOFTWARE_VERSIONS        } from './modules/goodwright/dump_software_versions/main'

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

    // Prepare manditory params into channels
    ch_samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_counts      = Channel.from(count_files)

    ch_meta = Channel.empty()
    if(params.run_input_check) {
        /*
        * MODULE: Check the samplesheet for errors
        */
        SAMPLE_DIFF_SAMPLESHEET_CHECK (
            ch_samplesheet,
            ch_counts.collect(),
            params.count_sep
        )
        ch_versions = ch_versions.mix(SAMPLE_DIFF_SAMPLESHEET_CHECK.out.versions)

        if (count_files.size() > 1) {
            ch_counts = SAMPLE_DIFF_SAMPLESHEET_CHECK.out.counts
        }
        // ch_counts | view

        /*
        * MODULE: Parse samplesheet into meta and fastq files
        */
        ch_meta = SAMPLE_DIFF_SAMPLESHEET_CHECK.out.csv
            .splitCsv ( header:true, sep:"," )
        //ch_meta | view
    }
    //EXAMPLE CHANNEL STRUCT: [sample_id:RAP1_IAA_30M_REP1, condition:B]
    //ch_meta | view
    //SAMPLE_DIFF_SAMPLESHEET_CHECK.out.csv | view

    ch_dsq_results = Channel.empty()
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
        *          but filter for only the conditions specified
        */
        ch_comparison_set = ch_meta
            .map { it[params.contrast_column] }
            .unique()
        //ch_comparison_set | view

        ch_comparisons = ch_comparison_set
            .combine(ch_comparison_set)
            .filter { it[0] != it[1] }
        //ch_comparisons | view

        if(comparisons && comparisons[0] != "all") {
            ch_comparisons = ch_comparisons
                .filter { ( it[0] + "_" + it[1] ) in comparisons }
        }
        //ch_comparisons | view

        ch_comparisons = ch_comparisons
                .map { [it[0] + "_" + it[1], it[0], it[1], '' ]}
        //ch_comparisons | view

        /*
        * CHANNEL: Join any set blocking factors with the target comparisons
        */
        if (ch_blocking_factors) {
            ch_comparisons = ch_comparisons
                .join (ch_blocking_factors, remainder: true)
                .map {
                    def block = ''
                    if(it[4]) {
                        block = it[4]
                    }
                    [it[0], it[1], it[2], block ] 
                }
        }
        //ch_comparisons | view

        /*
        * MODULE: Run deseq2
        */
        R_DESEQ2 (
            ch_design.collect(),
            ch_counts.collect(),
            params.contrast_column,
            ch_comparisons.map { it[1] },
            ch_comparisons.map { it[2] },
            ch_comparisons.map { it[3] }
        )
        ch_dsq_results = R_DESEQ2.out.results
        ch_versions = ch_versions.mix(R_DESEQ2.out.versions)
        //R_DESEQ2.out.results | view

        /*
        * CHANNEL: Get the first rdata object from deseq2
        */
        ch_dsq2_rdata = R_DESEQ2.out.rdata
            .collect()
            .map{[[id:'dsq2'], it[1]]}
        //ch_dsq2_rdata | view

        if (params.run_study_plots) {
            /*
            * MODULE: Run deseq2 plots
            */
            R_DESEQ2_PLOTS (
                ch_dsq2_rdata
            )
            ch_versions = ch_versions.mix(R_DESEQ2_PLOTS.out.versions)

            /*
            * MODULE: Run pcaexplorer
            */
            R_PCAEXPLORER (
                ch_dsq2_rdata
            )
            ch_versions = ch_versions.mix(R_PCAEXPLORER.out.versions)
        }

        /*
        * MODULE: Run Volcano Plot
        */
        if (params.run_volcano) {
            R_VOLCANO_PLOT (
                R_DESEQ2.out.results,
                params.contrast_column,
                ch_comparisons.map { it[1] },
                ch_comparisons.map { it[2] },
                params.blocking_factors
            )
            ch_versions = ch_versions.mix(R_VOLCANO_PLOT.out.versions)
        }
    }

    if(params.run_gsea) {
        /*
        * MODULE: Prep GSEA input channel with the filename as the meta id
        */
        ch_gsea_input = ch_dsq_results
            .map{ [[id: it[1].simpleName], it[1]] }
        //ch_gsea_input | view

        /*
        * MODULE: Run GSEA
        */
        R_GSEA (
            ch_gsea_input,
            params.organism
        )
        ch_versions = ch_versions.mix(R_GSEA.out.versions)
    }

    /*
    * MODULE: Collect software versions
    */
    DUMP_SOFTWARE_VERSIONS (
        ch_versions.unique().collectFile()
    )
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
