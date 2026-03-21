#!/usr/bin/env nextflow

process PICARD_CREATE_DICT {
    label 'process_single'
    conda 'envs/picard_env.yml'
    tag "PICARD_CREATE_DICT"
    publishDir "${params.outdir}/picard_testing", mode: 'copy'
    
    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict

    """
    picard CreateSequenceDictionary R=$fasta O=${fasta.baseName}.dict
    """
}