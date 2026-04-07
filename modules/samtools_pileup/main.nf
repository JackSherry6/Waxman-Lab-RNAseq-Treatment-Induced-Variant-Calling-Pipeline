
process SAMTOOLS_PILEUP {
    label 'process_single'
    conda 'envs/samtools_env.yml'
    //publishDir params.outdir, mode:'copy'  //mode:'move'
    tag "control_vs_${exp_name}" 

    input: 
    tuple val(exp_name), path(exp_bam), path(exp_bai), val(control_name), path(control_bam), path(control_bai)
    //tuple path(bam1), path(bai1), path(bam2), path(bai2)
    path(ref_genome)

    output:
    tuple val("control_vs_${exp_name}"), path("control_vs_${exp_name}.mpileup"), val(exp_name)

    script:
    """
    samtools mpileup -f $ref_genome -q 1 -B $control_bam $exp_bam > control_vs_${exp_name}.mpileup
    """
    // -q 1 is min mapping quality of 1. -B is to disable BAQ adjustment which is recommended for varscan

    stub:
    """
    touch control_vs_${exp_name}.mpileup
    """
}
