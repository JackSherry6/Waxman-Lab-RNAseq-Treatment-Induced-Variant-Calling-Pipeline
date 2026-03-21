process VARSCAN_TUMOR {
    label 'process_single'
    conda 'bioconda::varscan'  // For some reason scc doesn't like varscan as a yml file
    tag "${mpileup.baseName}"
    //publishDir params.outdir, mode:'copy'

    input: 
    tuple val(sample), path(mpileup), val(exp)
    val (control_cnt)

    output:
    tuple val("${sample}_snp"), path("${sample}_out*.snp.vcf"), emit : snp
    tuple val("${sample}_indel"), path("${sample}_out*.indel.vcf"), emit : indel

    script:
    """
    varscan somatic $mpileup ${sample}_out \
        --mpileup 1 \
        --min-coverage 10 \
        --min-coverage-normal $control_cnt \
        --min-coverage-tumor 10 \
        --min-var-freq 0.10 \
        --min-freq-for-hom 0.75 \
        --normal-purity 0.95 \
        --tumor-purity 0.8 \
        --somatic-p-value 0.05 \
        --strand-filter 1 \
        --output-vcf 1
    """

    stub:
    """
    touch ${sample}_out.snp.vcf
    touch ${sample}_out.indel.vcf
    """
}
/*
varscan somatic $mpileup ${sample}_out \
        --mpileup 1 \
        --min-coverage 10 \
        --min-var-freq 0.10 \
        --somatic-p-value 0.05 \
        --output-vcf 1
        */