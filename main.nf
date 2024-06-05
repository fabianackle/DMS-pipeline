#!/usr/bin/env nextflow

/* Define parameters */
params.reads = "$projectDir/data/*_R{1,2}_*"
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
    path "${sample_id}_R1_adapter_removed.fastq.gz"
    path "${sample_id}_R2_adapter_removed.fastq.gz"

    script:
    """
    cutadapt  -j 0 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o ${sample_id}_R1_adapter_removed.fastq.gz \
        -p ${sample_id}_R2_adapter_removed.fastq.gz \
        ${reads[0]} ${reads[1]}
    """
}

/* Workflow */
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    RemoveAdapter(read_pairs_ch)
}