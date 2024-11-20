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
    tuple val(sample_id), path("${sample_id}_R1_adapter_removed.fastq.gz"), path("${sample_id}_R2_adapter_removed.fastq.gz"), emit: cut
    path("${sample_id}_cutadapt.log"), emit: log

    script:
    """
    cutadapt -j $task.cpus \
        -a ${params.adapter_R1} \
        -A ${params.adapter_R2} \
        -o ${sample_id}_R1_adapter_removed.fastq.gz \
        -p ${sample_id}_R2_adapter_removed.fastq.gz \
        --minimum-length 50 \
        ${reads[0]} ${reads[1]} \
        > ${sample_id}_cutadapt.log
    """

    stub:
    """
    touch ${sample_id}_R1_adapter_removed.fastq.gz ${sample_id}_R2_adapter_removed.fastq.gz
    touch ${sample_id}_cutadapt.log
    """
}

process FastQC {
    cpus 8
    memory '4 GB'
    time '10m'
    conda "bioconda::fastqc=0.12.1"
    tag "FastQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.{html,zip}"), emit: stats

    script:
    """
    fastqc --threads $task.cpus \
        ${reads[0]} ${reads[1]}
    """

    stub:
    """
    touch ${sample_id}.html ${sample_id}.zip
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
    tuple val(sample_id), path("${sample_id}_aligned.raw.bam"), emit: aligned

    script:
    """
    bwa index $wt_sequence

    bwa mem -t $task.cpus \
        $wt_sequence \
        $trimmed_sequence_1 $trimmed_sequence_2 \
        | samtools view -S -b - > ${sample_id}_aligned.raw.bam
    """

    stub:
    """
    touch ${sample_id}_aligned.raw.bam
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
    path("${sample_id}_sorted.raw.bam"), emit: sorted

    script:
    """
    cp ${aligned_bam} ${sample_id}_sorted.raw.bam
    sambamba sort -t $task.cpus ${sample_id}_sorted.raw.bam
    """
    stub:
    """
    touch ${sample_id}_sorted.raw.bam
    """
}

process AlignSort {
    cpus 8
    memory '8 GB'
    time '45m'
    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    tag "BWA and samtools on $sample_id"

    input:
    tuple val(sample_id), path(trimmed_sequence_1), path(trimmed_sequence_2), path(wt_sequence)

    output:
    tuple val(sample_id), path("${sample_id}_adaptor_removed_trimmed.raw.bam"), emit: bam

    script:
    """
    bwa index $wt_sequence

    bwa mem -t $task.cpus \
        $wt_sequence \
        $trimmed_sequence_1 $trimmed_sequence_2 \
        | samtools sort -@ $task.cpus \
        -o ${sample_id}_adaptor_removed_trimmed.raw.bam
    """

    stub:
    """
    touch ${sample_id}_adaptor_removed_trimmed.raw.bam
    """   
}

process Subsample {
    cpus 8
    memory '1 GB'
    time '10m'
    conda "bioconda::samtools=1.20"
    tag "Samtools view on $sample_id"

    input:
    tuple val(sample_id), path(big_bam)

    output:
    tuple val(sample_id), path("${sample_id}_subsampled.bam"), emit: bam

    script:
    """
    samtools view -@ $task.cpus \
        --subsample 0.01 --subsample-seed 123 \
        -b -o "${sample_id}_temp.bam" \
        $big_bam

    samtools sort -@ $task.cpus \
        "${sample_id}_temp.bam" \
        -o "${sample_id}_subsampled.bam"

    rm "${sample_id}_temp.bam"
    """

    stub:
    """
    touch ${sample_id}_subsampled.bam
    """
}

process SamtoolsStats {
    cpus 8
    memory '1 GB'
    time '45m'
    conda "bioconda::samtools=1.20"
    tag "Samtools on $sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    path("*.txt"), emit: stats

    script:
    """
    samtools stats -@ $task.cpus \
        ${bam} > ${sample_id}_samtools_stats.txt
    samtools coverage \
        ${bam} > ${sample_id}_samtools_coverage.txt
    """

    stub:
    """
    touch ${sample_id}_samtools_stats.txt
    touch ${sample_id}_samtools_coverage.txt
    """
}

process SamtoolsStats2 {
    cpus 8
    memory '1 GB'
    time '45m'
    conda "bioconda::samtools=1.20"
    tag "Samtools on $sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    path("*.txt"), emit: stats

    script:
    """
    samtools stats -@ $task.cpus \
        ${bam} > ${sample_id}_samtools_stats_dms.txt
    samtools coverage \
        ${bam} > ${sample_id}_samtools_coverage_dms.txt
    """

    stub:
    """
    touch ${sample_id}_samtools_stats_dms.txt
    touch ${sample_id}_samtools_coverage_dms.txt
    """
}

process Analysis_DMS {
    cpus 1
    memory '1 GB'
    time '2h'
    conda "DMS_ABC.yml"
    tag "DMS_ABC on $sample_id"

    publishDir params.outdir, mode: 'copy', pattern: "*.{csv,txt}"

    input:
    tuple val(sample_id), path(bam), path(wt_sequence)

    output:
    path("${sample_id}_triplet_count.txt_readingframe_{1,2}_HDF5.csv"), emit: counts_readingframe
    path("${sample_id}_triplet_count.txt"), emit: counts
    tuple val(sample_id), path("${sample_id}_codontruncated.bam"), emit: bam
    
    script:
    """
    run_dms_abc.py \
        --sample_id ${sample_id} \
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

    stub:
    """
    touch ${sample_id}_triplet_count.txt_readingframe_1_HDF5.csv ${sample_id}_triplet_count.txt_readingframe_2_HDF5.csv
    touch ${sample_id}_triplet_count.txt
    touch ${sample_id}_codontruncated.bam
    """
}

process MultiQC {
    cpus 1
    time '2m'
    memory '256 MB'
    conda "bioconda::multiqc=1.25.1"
    tag "MultiQC"

    publishDir params.outdir, mode: 'copy'

    input:
    path("*")
    path(multiqc_config)

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data

    script:
    """
    multiqc --config ${multiqc_config} .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html
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

    if (workflow.stubRun) {
        // Stub channels
        read_pairs_ch = Channel.fromList([
            ['sample1', [file('dummy_R1.fastq.gz'), file('dummy_R2.fastq.gz')]],
            ['sample2', [file('dummy_R1.fastq.gz'), file('dummy_R2.fastq.gz')]]
        ])
        wt_sequence_ch = Channel.fromList([file('dummy_wt_sequence.fasta')])
    } else {
        // Real channels
        read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
        wt_sequence_ch = Channel.fromPath(params.wt_sequence, checkIfExists: true)
    }

    FastQC(read_pairs_ch)
    RemoveAdapter(read_pairs_ch)
    AlignSort(RemoveAdapter.out.cut.combine(wt_sequence_ch))
    SamtoolsStats(AlignSort.out.bam)
    Analysis_DMS(AlignSort.out.bam.combine(wt_sequence_ch))
    SamtoolsStats2(Analysis_DMS.out.bam)

    multiqc_config_ch = Channel.fromPath(params.multiqc_config)
    multiqc_files_ch = Channel.empty()
    multiqc_files_ch = multiqc_files_ch.mix(FastQC.out.stats)
    multiqc_files_ch = multiqc_files_ch.mix(RemoveAdapter.out.log)
    multiqc_files_ch = multiqc_files_ch.mix(SamtoolsStats.out.stats)
    multiqc_files_ch = multiqc_files_ch.mix(SamtoolsStats2.out.stats)
    multiqc_files_ch = multiqc_files_ch.collect()
    MultiQC(multiqc_files_ch, multiqc_config_ch)
}