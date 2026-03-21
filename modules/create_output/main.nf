process CREATE_OUTPUT {
    label 'process_single'
    conda 'envs/biopython_env.yml'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(group), path(indel_csvs), path(snp_csvs)

    output:
    path "${group}_snps.csv"
    path "${group}_indels.csv"

    script:
    """
    set -euo pipefail

    group_variants.py \\
      --group ${group} \\
      --indels ${indel_csvs.join(' ')} \\
      --snps   ${snp_csvs.join(' ')} \\
      -o .
    """
}

// groups each of the samples per group into a combined output csv file