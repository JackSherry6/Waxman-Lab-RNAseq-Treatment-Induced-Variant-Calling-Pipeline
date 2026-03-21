process PROCESS_VCFS {
    label 'process_single'
    conda 'envs/biopython_env.yml'
    publishDir "${params.outdir}/individual_csvs", mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path(lncrna_ref)
    path(known_snps)
    
    output:
    tuple val(sample), path("${sample}.processed.csv")

    script:
    def lncrna_arg = (lncrna_ref.name != 'NO_FILE_LNCRNA') ? "--lncrna ${lncrna_ref}" : ""
    def snps_arg   = (known_snps.name != 'NO_FILE_KNOWNSNPS') ? "--known_snps ${known_snps}" : ""

    """
    process_vcfs.py ${vcf} ${sample}.processed.csv ${lncrna_arg} ${snps_arg}
    """

    stub:
    """
    touch ${sample}.processed.csv
    """
}

// Expands VCF into a more manageable CSV format:
// - Splits INFO field into individual columns
// - Parses SnpEff ANN annotations (if present)
// - Splits FORMAT/NORMAL/TUMOR into sample-specific columns
// - Adds lncRNA column if variant overlaps lncRNA exon coordinates
// - Adds known_snp_* columns if variant position exists in consensus VCF
//
// Optional inputs:
// - lncrna_ref: GTF-like file with lncRNA exon coordinates
// - known_snps: Consensus VCF (e.g., from FreeBayes) for cohort variants
