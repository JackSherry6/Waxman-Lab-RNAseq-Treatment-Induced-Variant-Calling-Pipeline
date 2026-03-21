process BED_TO_INTERVAL_LIST {
    label 'process_single'
    conda 'envs/picard_env.yml'
    tag "BED_TO_INTERVAL_LIST"
    publishDir "${params.outdir}/picard_testing", mode: 'copy'

    input:
    path bed
    path dict

    output:
    path "rrna.interval_list", emit: intervals

    """
    picard BedToIntervalList I=$bed O=rrna.interval_list SD=$dict
    """
}