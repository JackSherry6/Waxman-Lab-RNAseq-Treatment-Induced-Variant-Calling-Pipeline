process BUILD_SNPEFF_DB {
    label 'process_low'
    conda 'envs/snpeff_env.yml'
    publishDir params.outdir, mode: 'copy'

    input:
    path ref_genome
    path gtf

    output:
    path("snpEff_data")

    script:
    """
    set -euo pipefail

    GENOME_NAME="cd1_snpeff_db"

    mkdir -p snpEff_data/data/\${GENOME_NAME}

    cp ${ref_genome} snpEff_data/data/\${GENOME_NAME}/sequences.fa
    cp ${gtf}        snpEff_data/data/\${GENOME_NAME}/genes.gtf

    printf "%s\n" "\${GENOME_NAME}.genome : \${GENOME_NAME}" > snpEff_data/snpEff.config

    snpEff build -gtf22 -v -noCheckCds -noCheckProtein -c snpEff_data/snpEff.config \${GENOME_NAME}
    """

    stub:
    """
    mkdir -p snpEff_data
    """
}