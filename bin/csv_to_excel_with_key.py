#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Convert CSV to Excel with a legend sheet.")
    parser.add_argument("input_csv", help="Input CSV file")
    parser.add_argument("-o", "--output", help="Output Excel file (.xlsx)", default="output.xlsx")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)

    individuals_legend = {
        "CHROM": "Chromosome",
        "POS": "Position on chromosome",
        "REF": "Reference Allele",
        "ALT": "Alternate Alelle",
        "SS": "Somatic Status (0=Ref, 1=Germline, 2=Somatic, 3=LOH, 5=Unexplained)",
        "SSC": "Single sample confidence (phred scale confidence score based on SPV, read depth, allele balance (ALT vs REF), and strand support)",
        "SPV": "Single position variant (chance of seeing this many ALT reads at a position assuming random sequencing errors)",
        "FREQ_normal": "Frequency of the alternate allele among the control samples",
        "FREQ_exp": "Frequency of the alternate allele among the experimental or tumor samples",
        "DP4_normal": "High quality read counts for control samples (REF_forward, REF_reverse, ALT_forward, ALT_reverse)",
        "Depth_normal": "Total high quality read counts for control samples",
        "DP4_exp": "High quality read counts for experimental/tumor samples (REF_forward, REF_reverse, ALT_forward, ALT_reverse)",
        "Depth_exp": "Total high quality read counts for experimental/tumor samples",
        "LOF": "Loss-of-function prediction (Gene | GeneID | num_transcripts | percent_affected)",
        "NMD": "Nonsense-mediated decay prediction (same subfield format as LOF)",
        "snpeff_effect": "Predicted consequence term(s) from SnpEff",
        "snpeff_impact": "SnpEff impact classification: HIGH / MODERATE / LOW / MODIFIER",
        "snpeff_gene_name": "RefSeq mRNA accession number",
        "snpeff_gene_id": "Gene stable ID (common gene name)",
        "snpeff_feature_type": "Feature affected (transcript or intergenic region)",
        "snpeff_feature_id": "Transcript ID",
        "snpeff_biotype": "Transcript biotype (protein-coding, pseudogene, or nan)",
        "snpeff_rank": "Which exon in the gene the SNP is in",
        "snpeff_hgvs_c": "HGVS coding DNA notation",
        "snpeff_hgvs_p": "HGVS protein notation for predicted amino-acid change",
        "snpeff_cdna_poslen": "cDNA position / length summary",
        "snpeff_cds_poslen": "CDS position / total CDS length",
        "snpeff_aa_poslen": "Amino acid position / total protein length",
        "snpeff_distance": "Distance to the feature for upstream/downstream/intergenic effects",
        "snpeff_errors_warnings_info": "SnpEff warnings/info flags",
        "lncRNA": "Overlapping lncRNAs at the SNP position",
        "known_snp_overlap": "Whether variant exists in the known SNPs VCF (Yes or No)",
        "known_snp_ref": "Reference allele recorded in the known SNPs VCF",
        "known_snp_alt": "Alternate allele(s) in the known SNPs VCF",
        "known_snp_qual": "QUAL field from the known SNPs VCF entry",
        "known_snp_af": "Allele frequency in the known SNPs VCF",
        "known_snp_dp": "Read depth in the known SNPs VCF",
        "known_snp_type": "Variant type (snp, ins, del)",
        "known_snp_ns": "Number of samples with data at that site",
        "known_snp_match": "Match type: exact / same_position_same_ref / same_position_diff_ref"
    }

    group_legend = {
    }

    import os
    basename = os.path.splitext(os.path.basename(args.input_csv))[0]
    if basename.endswith('_snps') or basename.endswith('_indels'):
        legend = group_legend
    else:
        legend = individuals_legend

    legend_df = pd.DataFrame(list(legend.items()), columns=["Column", "Description"])

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="data", index=False)
        legend_df.to_excel(writer, sheet_name="legend", index=False)

    print(f"Excel file written to: {args.output}")

if __name__ == "__main__":
    main()
