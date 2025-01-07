process SUBSAMPLE {
    cpus 8
    memory '1 GB'
    time '10m'
    conda "bioconda::samtools=1.20"
    tag "$sample_id"

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