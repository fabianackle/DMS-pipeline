process REMOVE_ADAPTER {
    errorStrategy 'ignore'
    cpus 8
    memory '1 GB'
    time '10m'
    conda "bioconda::cutadapt=4.8"
    tag "$sample_id"

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