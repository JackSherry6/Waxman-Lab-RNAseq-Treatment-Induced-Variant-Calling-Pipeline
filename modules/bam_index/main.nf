#!/usr/bin/env nextflow

process BAM_INDEX {
    label 'process_single'
    conda 'envs/samtools_env.yml'
    tag "${name}"
    publishDir "${params.outdir}/bams", mode: 'copy'

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path(bam), path("${bam}.bai")

    script:
    """
    samtools index $bam
    """

    stub:
    """
    touch ${bam}.bai
    """

}