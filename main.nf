#!/usr/bin/env nextflow

/* Define parameters */
params.reads = "$projectDir/data/*_R{1,2}*.fastq.gz"
params.wt_sequence = "$projectDir/data/EfrEF_opt_wt_sequence.fa"
params.outdir = "$projectDir/results"

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
    errorStrategy 'ignore'
    cpus 8
    conda "bioconda::cutadapt=4.8"
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
    cpus 8
    conda "bioconda::bwa-mem2=2.2.1"
    tag "BWA on $sample_id"

    input:
    tuple path(wt_sequence), val(sample_id), path(trimmed_sequence_1), path(trimmed_sequence_2)

    output:
    tuple val(sample_id), path("${sample_id}_adaptor_removed_trimmed.raw.sam")

    script:
    """
    bwa-mem2 index $wt_sequence

    bwa-mem2 mem -t $task.cpus \
        $wt_sequence \
        $trimmed_sequence_1 $trimmed_sequence_2 \
        > ${sample_id}_adaptor_removed_trimmed.raw.sam
    """
}

process Sort {
    cpus 8
    conda "bioconda::samtools=1.20"
    tag "Samtools on $sample_id"

    input:
    tuple val(sample_id), path(aligned_sam)

    output:
    path("${sample_id}_adaptor_removed_trimmed.raw.bam")

    script:
    """
    samtools sort -@ $task.cpus \
        $aligned_sam \
        -o ${sample_id}_adaptor_removed_trimmed.raw.bam
    """
}

process Subsample {
    cpus 8
    conda "bioconda::samtools=1.20"
    tag "Samtools view on $big_bam"

    publishDir params.outdir, mode: 'copy'

    input:
    path big_bam

    output:
    path "${big_bam.baseName}_subsampled.bam"

    script:
    """
    samtools view -@ $task.cpus \
        --subsample 0.01 --subsample-seed 123 \
        -b -o "${big_bam.baseName}_temp.bam" \
        $big_bam

    samtools sort -@ $task.cpus \
        "${big_bam.baseName}_temp.bam" \
        -o "${big_bam.baseName}_subsampled.bam"

    rm "${big_bam.baseName}_temp.bam"
    """
}

process Analysis_DMS {
    cpus 1
    conda "DMS_ABC.yml"
    tag "DMS_ABC on $bam"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple path(wt_sequence), path(bam)

    output:
    path "${bam.baseName}_triplet_count.txt_readingframe_{1,2}_HDF5.csv"
    path "${bam.baseName}_triplet_count.txt"
    

    script:
    """
    run_dms_abc.py \
        --bam "${bam}" \
        --reference "${wt_sequence}" \
        --positions 72 84 207 216 360 372 381 384 396 405 519 528 540 552 672 684 696 705 717 729 741 750 762 825 828 837 840 849 858 861 870 879 882 891 1026 1119 1245 1419 1488 1506 1581 1879 1891 2050 2059 2203 2218 2224 2227 2239 2248 2362 2371 2383 2395 2515 2527 2539 2545 2548 2560 2572 2584 2593 2605 2671 2680 2683 2692 2701 2704 2713 2722 2725 2734 2863 2953 3079 3253 3322 3340 3415 \
        --wt_ref_position 2545 \
        --wt_codon TTG \
        --wt_count 10000 \
        --readingframes \
        --frameshift_position 1728 \
        --frameshift_offset 52 \
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
    align_input_ch = wt_sequence_ch.combine(trimmed_ch)
    aligned_ch = Align(align_input_ch)
    sorted_ch = Sort(aligned_ch)
    //subsample_ch = Subsample(sorted_ch)
    //dms_abc_input_ch = wt_sequence_ch.combine(subsample_ch)
    dms_abc_input_ch = wt_sequence_ch.combine(sorted_ch)
    //dms_abc_input_ch.view()
    Analysis_DMS(dms_abc_input_ch)
}