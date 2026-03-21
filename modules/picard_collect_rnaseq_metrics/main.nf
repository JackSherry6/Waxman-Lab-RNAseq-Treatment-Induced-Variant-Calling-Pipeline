process PICARD_COLLECT_RNASEQ_METRICS {
    label 'process_single'
    conda 'envs/picard_env.yml'
    tag "PICARD_COLLECT_RNASEQ_METRICS"

    input:
    tuple val(sample), path(bam)
    path refflat
    path rrna_intervals
    path ref_fasta

    output:
    path "${sample}.picard.rna_metrics.txt", emit: metrics

    """
    picard CollectRnaSeqMetrics \
        I=$bam \
        O=${sample}.picard.rna_metrics.txt \
        REF_FLAT=$refflat \
        RIBOSOMAL_INTERVALS=$rrna_intervals \
        STRAND_SPECIFICITY=NONE \
        R=$ref_fasta
    """
}