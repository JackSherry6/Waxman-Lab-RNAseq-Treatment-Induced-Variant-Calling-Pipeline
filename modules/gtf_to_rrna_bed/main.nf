process GTF_TO_RRNA_BED {
    label 'process_single'
    tag "GTF_TO_RRNA_BED"
    publishDir "${params.outdir}/picard_testing", mode: 'copy'

    input:
    path gtf

    output:
    path "rrna.bed", emit: bed

    """
    awk 'BEGIN{OFS="\\t"}
        \$0 !~ /^#/ && \$3=="exon" && (\$0 ~ /rRNA/ || \$0 ~ /rrna/) {
            chr=\$1; start=\$4-1; end=\$5; strand=\$7;
            print chr,start,end,"rRNA",0,strand
        }' "${gtf}" | sort -k1,1 -k2,2n > rrna.bed
    """
}