process MERGE_BAMS {
    label 'process_low'
    conda 'envs/samtools_env.yml'
    tag "${sample}"
    //publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample), path(bams), path(bais)
    
    output:
    tuple val(sample), path("${sample}.merged.bam"), path("${sample}.merged.bam.bai")
    
    script:
    """
    samtools merge -@ $task.cpus ${sample}.merged.bam ${bams}
    samtools index -@ $task.cpus ${sample}.merged.bam
    """

    stub:
    """
    touch ${sample}.merged.bam
    touch ${sample}.merged.bam.bai
    """
}