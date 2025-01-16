process ALIGN_SORT {
    cpus 8
    memory '8 GB'
    time '1h'
    conda "bioconda::bwa=0.7.18 bioconda::samtools=1.20"
    tag "$sample_id"

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