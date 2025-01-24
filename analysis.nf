#!/usr/bin/env nextflow

/* Processes */
include { INGEST                             } from './modules/ingest.nf'
include { ANALYSIS_DMS                       } from './modules/analysis_dms.nf'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_1 } from './modules/samtools_stats.nf'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_2 } from './modules/samtools_stats.nf'
include { MULTIQC                            } from './modules/multiqc.nf'

/* Workflows */
workflow dms {
    take:
    bams_ch
    wt_sequence_ch 
    multiqc_config_ch 

    main:
    multiqc_files_ch = Channel.empty()

    INGEST(bams_ch)

    SAMTOOLS_STATS_1(INGEST.out.bam, 'aligned')
    multiqc_files_ch = multiqc_files_ch.mix(SAMTOOLS_STATS_1.out.stats)

    ANALYSIS_DMS(INGEST.out.bam.combine(wt_sequence_ch))
    SAMTOOLS_STATS_2(ANALYSIS_DMS.out.bam, 'dms')
    multiqc_files_ch = multiqc_files_ch.mix(SAMTOOLS_STATS_2.out.stats)

    multiqc_files_ch = multiqc_files_ch.collect()
    MULTIQC(multiqc_files_ch, multiqc_config_ch)    
}

workflow {
    /* Print pipeline info */
    log.info """
    ===============================
    D M S - A B C   P I P E L I N E
    ===============================
    Output dir: ${params.outdir}
    """
    .stripIndent()

    if (workflow.stubRun) {
        // Stub channels
        bams_ch = Channel.fromList([
            file('file0.bam'),
            file('file2.bam'),
            file('file3.bam'),
            ])
        wt_sequence_ch = Channel.fromList([file('dummy_wt_sequence.fasta')])
    } else {
        // Real channels
        bams_ch = Channel.fromFilePairs(params.bams, checkIfExists: true)
        wt_sequence_ch = Channel.fromPath(params.wt_sequence, checkIfExists: true)
    }

    multiqc_config_ch = Channel.fromPath(params.multiqc_config)

    dms(bams_ch, wt_sequence_ch, multiqc_config_ch)
}