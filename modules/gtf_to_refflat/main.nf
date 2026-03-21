
process GTF_TO_REFFLAT {
    label 'process_single'
    tag "GTF_TO_REFFLAT"
    conda "envs/refflat_env.yml"
    publishDir "${params.outdir}/picard_testing", mode: 'copy'

    input:
    path gtf

    output:
    path "genes.refFlat", emit: refflat

    script:
    """
    gtfToGenePred "${gtf}" genes.genePred

    # refFlat = genePred with an extra first column (geneName). Here we use transcript name for geneName.
    awk 'BEGIN{OFS="\\t"} {print \$1,\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' genes.genePred > genes.refFlat
    """
}