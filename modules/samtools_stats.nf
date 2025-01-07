process SAMTOOLS_STATS {
    cpus 8
    memory '1 GB'
    time '45m'
    conda "bioconda::samtools=1.20"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)
    val(stage)

    output:
    path("*.txt"), emit: stats

    script:
    """
    samtools stats -@ $task.cpus \
        ${bam} > ${sample_id}_samtools_stats_${stage}.txt
    samtools coverage \
        ${bam} > ${sample_id}_samtools_coverage_${stage}.txt
    """

    stub:
    """
    touch ${sample_id}_samtools_stat_${stage}s.txt
    touch ${sample_id}_samtools_coverage_${stage}.txt
    """
}