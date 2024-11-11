#!/usr/bin/env nextflow

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
        -a ${params.adapter_R1} \
        -A ${params.adapter_R2} \
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
    time '1h'
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
    time '1h'
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
        --positions ${params.positions} \
        --nnk_positions ${params.nnk_positions} \
        --wt_ref_position ${params.wt_ref_position} \
        --wt_codon ${params.wt_codon} \
        --wt_count ${params.wt_count} \
        --readingframes \
        --frameshift_position ${params.frameshift_position} \
        --frameshift_offset ${params.frameshift_offset}
    """
}

/* Workflow */
workflow {
    /* Print pipeline info */
    log.info """
    ===============================
    D M S - A B C   P I P E L I N E
    ===============================
    Reads : ${params.reads}
    Reference: ${params.wt_sequence}
    Output dir: ${params.outdir}
    """
    .stripIndent()

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