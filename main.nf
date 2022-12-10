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

include { summary_log     } from './modules/goodwright/util/logging/main'
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

// // Check manditory input parameters to see if the files exist if they have been specified
// check_param_list = [
//     samplesheet: params.samplesheet,
//     fasta: params.fasta,
//     smrna_fasta: params.smrna_fasta,
//     gtf: params.gtf
// ]
// for (param in check_param_list) { 
//     if (!param.value) { 
//         exit 1, "Required parameter not specified: ${param.key}"
//     } 
//     else {
//         file(param.value, checkIfExists: true)
//     }
// }

// // Check non-manditory input parameters to see if the files exist if they have been specified
// check_param_list = [
//     params.target_genome_index,
//     params.smrna_genome_index
// ]
// for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

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

    // // Prepare manditory params into channels 
    // ch_samplesheet = file(params.samplesheet, checkIfExists: true)
    // ch_fasta       = file(params.fasta, checkIfExists: true)
    // ch_smrna_fasta = file(params.smrna_fasta, checkIfExists: true)
    // ch_gtf         = file(params.gtf, checkIfExists: true)

    // // Prepare non-manditory params into channels
    // if(params.target_genome_index) { ch_target_genome_index = file(params.target_genome_index, checkIfExists: true) }
    // if(params.smrna_genome_index)  { ch_smrna_genome_index = file(params.smrna_genome_index, checkIfExists: true) }

    // // Prepare genome and build indexes if required
    // ch_fasta_fai              = Channel.empty()
    // ch_filtered_gtf           = Channel.empty()
    // ch_chrom_sizes            = Channel.empty()
    // ch_smrna_fasta_fai        = Channel.empty()
    // ch_smrna_chrom_sizes      = Channel.empty()
    // ch_longest_transcript     = Channel.empty()
    // ch_seg_gtf                = Channel.empty()
    // ch_seg_filt_gtf           = Channel.empty()
    // ch_seg_resolved_gtf       = Channel.empty()
    // ch_seg_resolved_gtf_genic = Channel.empty()
    // if (params.run_genome_prep) {
    //     /*
    //     * SUBWORKFLOW: Prepare clipseq genome files
    //     */
    //     PREPARE_CLIPSEQ (
    //         ch_fasta,
    //         ch_smrna_fasta,
    //         ch_gtf,
    //         ch_target_genome_index,
    //         ch_smrna_genome_index
    //     )
    //     ch_versions               = ch_versions.mix(PREPARE_CLIPSEQ.out.versions)
    //     ch_fasta                  = PREPARE_CLIPSEQ.out.fasta
    //     ch_fasta_fai              = PREPARE_CLIPSEQ.out.fasta_fai
    //     ch_gtf                    = PREPARE_CLIPSEQ.out.gtf
    //     ch_filtered_gtf           = PREPARE_CLIPSEQ.out.filtered_gtf
    //     ch_chrom_sizes            = PREPARE_CLIPSEQ.out.chrom_sizes
    //     ch_smrna_fasta            = PREPARE_CLIPSEQ.out.smrna_fasta
    //     ch_smrna_fasta_fai        = PREPARE_CLIPSEQ.out.smrna_fasta_fai
    //     ch_smrna_chrom_sizes      = PREPARE_CLIPSEQ.out.smrna_chrom_sizes
    //     ch_longest_transcript     = PREPARE_CLIPSEQ.out.longest_transcript
    //     ch_seg_gtf                = PREPARE_CLIPSEQ.out.seg_gtf
    //     ch_seg_filt_gtf           = PREPARE_CLIPSEQ.out.seg_filt_gtf
    //     ch_seg_resolved_gtf       = PREPARE_CLIPSEQ.out.seg_resolved_gtf
    //     ch_seg_resolved_gtf_genic = PREPARE_CLIPSEQ.out.seg_resolved_gtf_genic
    //     ch_target_genome_index    = PREPARE_CLIPSEQ.out.genome_index
    //     ch_smrna_genome_index     = PREPARE_CLIPSEQ.out.smrna_index
    // }

    // ch_fastq = Channel.empty()
    // if(params.run_input_check) {
    //     /*
    //     * SUBWORKFLOW: Read in samplesheet, validate, stage input files and merge replicates
    //     */
    //     PARSE_FASTQ_INPUT (
    //         ch_samplesheet
    //     )
    //     ch_versions = ch_versions.mix(PARSE_FASTQ_INPUT.out.versions)
    //     ch_fastq    = PARSE_FASTQ_INPUT.out.fastq
    // }
    // //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false], [FASTQ]]
    // //ch_fastq | view



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