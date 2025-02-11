process SUBSAMPLE {
    conda "bioconda::samtools=1.20"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(big_bam)

    output:
    tuple val(sample_id), path("${sample_id}_subsampled.bam"), emit: bam

    script:
    """
    samtools view -@ $task.cpus \
        --subsample 0.01 \
        --subsample-seed 123 \
        -b $big_bam | \
    samtools sort -@ $task.cpus \
        -O bam \
        -o "${sample_id}_subsampled.bam"
    """

    stub:
    """
    touch ${sample_id}_subsampled.bam
    """
}