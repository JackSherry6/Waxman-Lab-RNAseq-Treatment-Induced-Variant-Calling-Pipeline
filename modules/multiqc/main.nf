process MULTIQC {
    label 'process_low'
    conda 'envs/multiqc_env.yml'
    publishDir "${params.outdir}/multiqc", mode: "copy"

    input:
    path(qc_files)

    output:
    path('*.html')
    path("multiqc_data")

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc ${qc_files} -f -o .
    """
}