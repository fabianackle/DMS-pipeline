#!/usr/bin/env nextflow

/* Define parameters */
params.reads = "$projectDir/data/*_R{1,2}_001.fastq.gz"
params.wt_sequence = "$projectDir/data/LmrCD_WT.fa"
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
    memory '1 GB'
    time '10m'
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
        --minimum-length 50 \
        ${reads[0]} ${reads[1]}
    """
}

process Align {
    cpus 8
    memory '2 GB'
    time '20m'
    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    tag "BWA on $sample_id"

    input:
    tuple path(wt_sequence), val(sample_id), path(trimmed_sequence_1), path(trimmed_sequence_2)

    output:
    tuple val(sample_id), path("${sample_id}_aligned.raw.bam")

    script:
    """
    bwa index $wt_sequence

    bwa mem -t $task.cpus \
        $wt_sequence \
        $trimmed_sequence_1 $trimmed_sequence_2 \
        | samtools view -S -b - > ${sample_id}_aligned.raw.bam
    """
}

process Sort {
    cpus 8
    memory '3 GB'
    time '10m'
    conda "bioconda::sambamba=1.0.1"
    tag "Sambamba on $sample_id"

    input:
    tuple val(sample_id), path(aligned_bam)

    output:
    path("${sample_id}_sorted.raw.bam")

    script:
    """
    cp ${aligned_bam} ${sample_id}_sorted.raw.bam
    sambamba sort -t $task.cpus ${sample_id}_sorted.raw.bam
    """
}

process AlignSort {
    cpus 8
    memory '8 GB'
    time '20m'
    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    tag "BWA and samtools on $sample_id"

    input:
    tuple path(wt_sequence), val(sample_id), path(trimmed_sequence_1), path(trimmed_sequence_2)

    output:
    path("${sample_id}_adaptor_removed_trimmed.raw.bam")

    script:
    """
    bwa index $wt_sequence

    bwa mem -t $task.cpus \
        $wt_sequence \
        $trimmed_sequence_1 $trimmed_sequence_2 \
        | samtools sort -@ $task.cpus \
        -o ${sample_id}_adaptor_removed_trimmed.raw.bam
    """    
}

process Subsample {
    cpus 8
    memory '1 GB'
    time '5 min'
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
    memory '1 GB'
    time '2h'
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
        --positions 66 78 216 225 369 381 390 393 405 414 528 537 549 561 681 693 705 714 726 738 750 759 771 846 849 858 861 870 879 882 891 900 903 912 1038 1131 1257 1431 1500 1518 1593 1863 1875 2217 2226 2259 2352 2370 2385 2391 2394 2406 2415 2529 2538 2550 2562 2580 2682 2694 2706 2715 2727 2739 2751 2760 2772 2853 2862 2865 2874 2883 2886 2895 2904 2907 2916 3042 3132 3258 3432 3501 3519 3594 \
        --nnk_positions 55 59 70 78 \
        --wt_ref_position 2259 \
        --wt_codon TCA \
        --wt_count 10000 \
        --readingframes \
        --frameshift_position 1740 \
        --frameshift_offset 3
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
    //aligned_ch = Align(align_input_ch)
    //sorted_ch = Sort(aligned_ch)
    sorted_ch = AlignSort(align_input_ch)
    subsample_ch = Subsample(sorted_ch)
    dms_abc_input_ch = wt_sequence_ch.combine(subsample_ch)
    //dms_abc_input_ch.view()
    Analysis_DMS(dms_abc_input_ch)
}