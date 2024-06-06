#!/usr/bin/env nextflow

/* Define parameters */
params.reads = "$projectDir/data/*_R{1,2}_*"
params.wt_sequence = "$projectDir/data/EfrEF_opt_wt_sequence.fa"
params.outdir = "results"

/* Print pipeline info */
log.info """
    ===============================
    D M S - A B C   P I P E L I N E
    ===============================
    Reads : ${params.reads}
    Output dir: ${params.outdir}
    """
    .stripIndent()

/* Processes */
process RemoveAdapter {
    tag "Cutadapt on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1_adapter_removed.fastq.gz"), path("${sample_id}_R2_adapter_removed.fastq.gz")

    script:
    """
    cutadapt -j $task.cpus \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o ${sample_id}_R1_adapter_removed.fastq.gz \
        -p ${sample_id}_R2_adapter_removed.fastq.gz \
        ${reads[0]} ${reads[1]}
    """
}

process Align {
    tag "BWA on $sample_id"

    input:
    path(wt_sequence)
    tuple val(sample_id), path(trimmed_sequence_1), path(trimmed_sequence_2)

    output:
    path "${sample_id}_adaptor_removed_trimmed.raw.bam"

    script:
    """
    bwa index $wt_sequence

    bwa mem -t $task.cpus \
        $wt_sequence \
        $trimmed_sequence_1 $trimmed_sequence_2 \
        > ${sample_id}_adaptor_removed_trimmed.raw.sam

    samtools sort -@ $task.cpus \
        ${sample_id}_adaptor_removed_trimmed.raw.sam \
        -o ${sample_id}_adaptor_removed_trimmed.raw.bam

    rm ${sample_id}_adaptor_removed_trimmed.raw.sam
    """
}

/* Workflow */
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    Channel
        .fromPath(params.wt_sequence, checkIfExists: true)
        .set { wt_sequence_ch }
    trimmed_ch = RemoveAdapter(read_pairs_ch)
    align_ch = Align(wt_sequence_ch, trimmed_ch)
}