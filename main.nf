#!/usr/bin/env nextflow

/* Processes */
include { REMOVE_ADAPTER                     } from './modules/remove_adapter.nf'
include { FASTQC                             } from './modules/fastqc.nf'
include { ALIGN_SORT                         } from './modules/align_sort.nf'
include { SUBSAMPLE                          } from './modules/subsample.nf'
include { ANALYSIS_DMS                       } from './modules/analysis_dms.nf'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_1 } from './modules/samtools_stats.nf'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_2 } from './modules/samtools_stats.nf'
include { MULTIQC                            } from './modules/multiqc.nf'

workflow {
    /* Print pipeline info */
    log.info(
        """
    ===============================
    D M S - A B C   P I P E L I N E
    ===============================
    Reads : ${params.reads}
    Reference: ${params.wt_sequence}
    Output dir: ${params.outdir}
    """.stripIndent()
    )

    if (workflow.stubRun) {
        // Stub channels
        read_pairs_ch = channel.fromList(
            [
                ['sample1', [file('dummy_R1.fastq.gz'), file('dummy_R2.fastq.gz')]],
                ['sample2', [file('dummy_R1.fastq.gz'), file('dummy_R2.fastq.gz')]],
            ]
        )
        wt_sequence_ch = channel.fromList([file('dummy_wt_sequence.fasta')])
    }
    else {
        // Real channels
        read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
        wt_sequence_ch = channel.fromPath(params.wt_sequence, checkIfExists: true)
    }

    multiqc_config_ch = channel.fromPath(params.multiqc_config)

    dms(read_pairs_ch, wt_sequence_ch, multiqc_config_ch)
}

/* Workflows */
workflow dms {
    take:
    read_pairs_ch
    wt_sequence_ch
    multiqc_config_ch

    main:
    multiqc_files_ch = channel.empty()
    FASTQC(read_pairs_ch)
    multiqc_files_ch = multiqc_files_ch.mix(FASTQC.out.stats)

    REMOVE_ADAPTER(read_pairs_ch)
    multiqc_files_ch = multiqc_files_ch.mix(REMOVE_ADAPTER.out.log)

    ALIGN_SORT(REMOVE_ADAPTER.out.cut.combine(wt_sequence_ch))
    SAMTOOLS_STATS_1(ALIGN_SORT.out.bam, 'aligned')
    multiqc_files_ch = multiqc_files_ch.mix(SAMTOOLS_STATS_1.out.stats)

    if (params.subsample) {
        SUBSAMPLE(ALIGN_SORT.out.bam)
        ANALYSIS_DMS(SUBSAMPLE.out.bam.combine(wt_sequence_ch))
    }
    else {
        ANALYSIS_DMS(ALIGN_SORT.out.bam.combine(wt_sequence_ch))
    }

    SAMTOOLS_STATS_2(ANALYSIS_DMS.out.bam, 'dms')
    multiqc_files_ch = multiqc_files_ch.mix(SAMTOOLS_STATS_2.out.stats)

    multiqc_files_ch = multiqc_files_ch.collect()
    MULTIQC(multiqc_files_ch, multiqc_config_ch)
}
