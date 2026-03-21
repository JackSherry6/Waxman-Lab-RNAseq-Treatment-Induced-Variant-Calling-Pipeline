process ANNOTATE_VARIANTS {
    label 'process_low'
    conda 'envs/snpeff_env.yml'
    publishDir "${params.outdir}/vcfs", mode: 'copy'

    input:
    tuple val(name), path(vcf)
    path(snpeff_db_dir)

    output:
    tuple val(name), path("${name}.snpeff.annotated.vcf")

    script:
    """
    set -euo pipefail

    GENOME_NAME="cd1_snpeff_db"

    snpEff -v \\
      -c ${snpeff_db_dir}/snpEff.config \\
      \${GENOME_NAME} \\
      ${vcf} > ${name}.snpeff.annotated.vcf
    """

    stub:
    """
    touch ${name}.snpeff.annotated.vcf
    """
}