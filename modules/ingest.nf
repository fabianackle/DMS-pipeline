process INGEST {
    cpus 8
    memory '1 GB'
    time '1m'
    tag "${bam.baseName}"

    input:
    path bam

    output:
    tuple val(bam.baseName), path(bam), emit: bam

    script:
    """
    """

    stub:
    """
    """
}