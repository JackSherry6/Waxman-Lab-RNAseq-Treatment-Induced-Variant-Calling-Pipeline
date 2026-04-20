process VARSCAN_CHR {
    label 'process_single'
    // Direct spec used because SCC does not resolve varscan yml files correctly
    conda 'bioconda::samtools=1.22.1 bioconda::varscan'
    tag "${exp_name}:${chrom}"

    input:
    tuple val(exp_name), path(exp_bam), path(exp_bai),
          val(ctrl_name), path(ctrl_bam), path(ctrl_bai),
          val(tumor_purity), val(chrom)
    path ref_genome
    val control_cnt

    output:
    tuple val(exp_name), path("${exp_name}_${chrom}.snp.vcf.gz"),   path("${exp_name}_${chrom}.snp.vcf.gz.tbi"),   emit: snp
    tuple val(exp_name), path("${exp_name}_${chrom}.indel.vcf.gz"), path("${exp_name}_${chrom}.indel.vcf.gz.tbi"), emit: indel

    script:
    """
    samtools mpileup -f ${ref_genome} -q 1 -B -r ${chrom} ${ctrl_bam} ${exp_bam} \
        > ${chrom}.mpileup

    varscan somatic ${chrom}.mpileup ${exp_name}_${chrom} \
        --mpileup 1 \
        --min-coverage 10 \
        --min-coverage-normal ${control_cnt} \
        --min-coverage-tumor 10 \
        --min-var-freq 0.10 \
        --min-freq-for-hom 0.75 \
        --normal-purity 0.95 \
        --tumor-purity ${tumor_purity} \
        --somatic-p-value 0.05 \
        --strand-filter 1 \
        --output-vcf 1

    rm -f ${chrom}.mpileup

    bgzip ${exp_name}_${chrom}.snp.vcf
    tabix -p vcf ${exp_name}_${chrom}.snp.vcf.gz
    bgzip ${exp_name}_${chrom}.indel.vcf
    tabix -p vcf ${exp_name}_${chrom}.indel.vcf.gz
    """

    stub:
    """
    touch ${exp_name}_${chrom}.snp.vcf.gz
    touch ${exp_name}_${chrom}.snp.vcf.gz.tbi
    touch ${exp_name}_${chrom}.indel.vcf.gz
    touch ${exp_name}_${chrom}.indel.vcf.gz.tbi
    """
}
