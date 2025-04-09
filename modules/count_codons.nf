process COUNT_CODONS {
    conda "bioconda::dnaio=1.2.2 bioconda::pysam=0.23.0 bioconda::samtools=1.21 conda-forge::polars=1.26"
    tag "$sample_id"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(wt_sequence)

    output:
    path("${sample_id}.csv"), emit: counts

    script:
    """
    count_codons.py \
        --sample_id ${sample_id} \
        --bam ${bam} \
        --reference ${wt_sequence} \
        --positions ${params.positions}
    """

    stub:
    """
    touch ${sample_id}.csv
    """
}