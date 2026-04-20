process MERGE_VCFS {
    label 'process_single'
    conda 'envs/bcftools_env.yml'
    tag "${sample}_${vcf_type}"

    input:
    tuple val(sample), val(vcf_type), path(vcfs), path(tbis)

    output:
    tuple val("${sample}_${vcf_type}"), path("${sample}_${vcf_type}.vcf")

    script:
    """
    bcftools concat -a --output-type v -o ${sample}_${vcf_type}.vcf ${vcfs}
    """

    stub:
    """
    touch ${sample}_${vcf_type}.vcf
    """
}
