process ANALYSIS_DMS {
    conda "DMS_ABC.yml"
    tag "$sample_id"

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