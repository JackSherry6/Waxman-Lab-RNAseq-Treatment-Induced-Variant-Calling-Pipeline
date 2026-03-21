process FILTER_VARIANTS {
    label 'process_single'
    conda 'envs/bcftools_env.yml'
    //publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample), path(vcf)
    
    output:
    tuple val(sample), path("${sample}.filtered.vcf")

    script:
    """
    bcftools view ${vcf} | \
    bcftools filter -i '(INFO/SS == 2 || INFO/SS == 3)' \
        -o ${sample}.filtered.vcf
    """

    stub:
    """
    touch ${sample}.filtered.vcf
    """
}
