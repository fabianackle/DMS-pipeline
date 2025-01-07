process FASTQC {
    cpus 8
    memory '4 GB'
    time '10m'
    conda "bioconda::fastqc=0.12.1"
    tag "$sample_id"

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