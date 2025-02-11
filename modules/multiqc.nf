process MULTIQC {
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