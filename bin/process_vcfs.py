#!/usr/bin/env python3
"""
VCF Processing Script
- Parses VCF files and expands INFO/FORMAT fields
- Annotates with lncRNA overlaps (if reference provided)
- Annotates with consensus/known SNP VCF overlaps (if provided)
"""

import re
import gzip
import argparse
import pandas as pd
from collections import defaultdict


def parse_snpeff_ann_first(ann_value: str) -> pd.Series:
    """Parse the first (most impactful) SnpEff annotation."""
    out = {
        "snpeff_allele": None, "snpeff_effect": None, "snpeff_impact": None,
        "snpeff_gene_name": None, "snpeff_gene_id": None, "snpeff_feature_type": None,
        "snpeff_feature_id": None, "snpeff_biotype": None, "snpeff_rank": None,
        "snpeff_hgvs_c": None, "snpeff_hgvs_p": None, "snpeff_cdna_poslen": None,
        "snpeff_cds_poslen": None, "snpeff_aa_poslen": None, "snpeff_distance": None,
        "snpeff_errors_warnings_info": None,
    }

    if pd.isna(ann_value):
        return pd.Series(out)

    ann_str = str(ann_value).strip()
    if not ann_str or ann_str == ".":
        return pd.Series(out)

    fields = ann_str.split(",")[0].split("|")
    fields += [""] * (17 - len(fields))

    keys = list(out.keys())
    for i, key in enumerate(keys):
        out[key] = fields[i] or None

    return pd.Series(out)


def read_vcf_to_df(vcf_path: str) -> pd.DataFrame:
    """Read a VCF file into a pandas DataFrame."""
    skip_rows = 0
    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                skip_rows += 1
            else:
                break

    df = pd.read_csv(vcf_path, sep="\t", skiprows=skip_rows, compression="infer")
    if "#CHROM" in df.columns:
        df = df.rename(columns={"#CHROM": "CHROM"})
    return df


def split_info_column(df: pd.DataFrame) -> pd.DataFrame:
    """Split the INFO column into separate columns."""
    info_df = (
        df["INFO"]
        .astype(str)
        .str.split(";")
        .apply(lambda items: {
            k: (v[0] if v else True)
            for item in items if item
            for k, *v in [item.split("=", 1)]
        })
        .apply(pd.Series)
    )
    return pd.concat([df.drop(columns=["INFO"]), info_df], axis=1)


def parse_format_row(format_str: str, normal_str: str, tumor_str: str) -> pd.Series:
    """Parse FORMAT field and split NORMAL/TUMOR sample data."""
    keys = str(format_str).split(":")
    row = {}
    for k, n, t in zip(keys, str(normal_str).split(":"), str(tumor_str).split(":")):
        row[f"{k}_normal"] = n
        row[f"{k}_tumor"] = t
    return pd.Series(row)


# =============================================================================
# lncRNA Reference Loading and Overlap Checking
# =============================================================================

def parse_gtf_attributes(attr_string: str) -> dict:
    """Parse GTF attribute string into a dictionary."""
    attrs = {}
    # Match patterns like key "value" or key 'value'
    pattern = r'(\w+)\s*["\']([^"\']*)["\']'
    matches = re.findall(pattern, attr_string)
    for key, value in matches:
        attrs[key] = value
    return attrs


def load_lncrna_reference(lncrna_file: str) -> dict:
    """
    Load lncRNA reference file (GTF-like format with exon coordinates).
    Returns a dictionary organized by chromosome for efficient lookup.
    
    Expected format (tab-separated):
    chr1  source  exon  start  end  score  strand  frame  attributes
    
    Attributes should contain gene_name or gene_id.
    """
    if not lncrna_file:
        return None
    
    print(f"Loading lncRNA reference from {lncrna_file}...")
    
    lncrna_regions = defaultdict(list)
    
    try:
        with open(lncrna_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                
                chrom = parts[0]
                # parts[1] = source, parts[2] = feature (exon)
                start = int(parts[3])
                end = int(parts[4])
                # parts[5] = score, parts[6] = strand, parts[7] = frame
                attributes = parts[8]
                
                # Parse attributes to get gene name
                attrs = parse_gtf_attributes(attributes)
                gene_name = attrs.get('gene_name', attrs.get('gene_id', 'unknown'))
                
                lncrna_regions[chrom].append({
                    'start': start,
                    'end': end,
                    'gene_name': gene_name
                })
        
        # Sort regions by start position for each chromosome
        for chrom in lncrna_regions:
            lncrna_regions[chrom].sort(key=lambda x: x['start'])
        
        total_regions = sum(len(v) for v in lncrna_regions.values())
        print(f"  Loaded {total_regions} lncRNA exon regions across {len(lncrna_regions)} chromosomes")
        
        return dict(lncrna_regions)
    
    except FileNotFoundError:
        print(f"  Warning: lncRNA file not found: {lncrna_file}")
        return None
    except Exception as e:
        print(f"  Warning: Error loading lncRNA file: {e}")
        return None


def check_lncrna_overlap(chrom: str, pos: int, lncrna_regions: dict) -> str:
    """
    Check if a SNP position overlaps with any lncRNA exon.
    Returns the gene name(s) if overlap found, empty string otherwise.
    """
    if not lncrna_regions or chrom not in lncrna_regions:
        return ""
    
    overlapping_genes = []
    
    for region in lncrna_regions[chrom]:
        # Check if position falls within the exon
        if region['start'] <= pos <= region['end']:
            overlapping_genes.append(region['gene_name'])
        # Early exit if we've passed the position (regions are sorted)
        elif region['start'] > pos:
            break
    
    # Return unique gene names, semicolon-separated
    unique_genes = list(dict.fromkeys(overlapping_genes))
    return ';'.join(unique_genes) if unique_genes else ""


def annotate_lncrna(df: pd.DataFrame, lncrna_regions: dict) -> pd.DataFrame:
    """Add lncRNA column to dataframe based on coordinate overlap."""
    if lncrna_regions is None:
        return df
    
    print("Annotating lncRNA overlaps...")
    
    df['lncRNA'] = df.apply(
        lambda row: check_lncrna_overlap(
            str(row['CHROM']),
            int(row['POS']),
            lncrna_regions
        ),
        axis=1
    )
    
    overlap_count = (df['lncRNA'] != '').sum()
    print(f"  Found {overlap_count} variants overlapping lncRNA exons")
    
    return df


# =============================================================================
# Consensus/Known SNPs VCF Loading and Overlap Checking
# =============================================================================

def load_known_snps_vcf(vcf_file: str) -> dict:
    """
    Load consensus/known SNPs VCF file and organize by (chrom, pos) for lookup.
    Returns a dictionary with (chrom, pos) as keys.
    
    Handles both FreeBayes and VarScan2 VCF formats.
    """
    if not vcf_file:
        return None
    
    print(f"Loading known SNPs VCF from {vcf_file}...")
    
    vcf_variants = {}
    sample_names = []
    
    try:
        # Handle gzipped or plain VCF
        opener = gzip.open if vcf_file.endswith('.gz') else open
        
        with opener(vcf_file, 'rt') as f:
            for line in f:
                line = line.strip()
                
                if not line:
                    continue
                
                # Parse header to get sample names
                if line.startswith('#CHROM'):
                    parts = line.split('\t')
                    if len(parts) > 9:
                        sample_names = parts[9:]
                    continue
                
                # Skip other header lines
                if line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 8:
                    continue
                
                chrom = parts[0]
                pos = int(parts[1])
                variant_id = parts[2]
                ref = parts[3]
                alt = parts[4]
                qual = parts[5]
                filter_val = parts[6]
                info = parts[7]
                
                # Parse INFO field into dict
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[item] = 'True'
                
                # Count how many samples have the variant (non ./. genotypes)
                samples_with_variant = 0
                sample_genotypes = []
                if len(parts) > 9:
                    for sample_data in parts[9:]:
                        gt = sample_data.split(':')[0] if sample_data else '.'
                        if gt not in ['.', './.', '.|.']:
                            sample_genotypes.append(gt)
                            # Count samples with at least one alt allele
                            if '1' in gt or '2' in gt:
                                samples_with_variant += 1
                
                key = (chrom, pos)
                vcf_variants[key] = {
                    'id': variant_id,
                    'ref': ref,
                    'alt': alt,
                    'qual': qual,
                    'filter': filter_val,
                    'af': info_dict.get('AF', ''),
                    'dp': info_dict.get('DP', ''),
                    'type': info_dict.get('TYPE', ''),
                    'ns': info_dict.get('NS', ''),           # Number of samples with data
                    'ac': info_dict.get('AC', ''),           # Allele count
                    'an': info_dict.get('AN', ''),           # Total alleles
                    'samples_with_variant': samples_with_variant,
                    'total_samples': len(sample_names)
                }
        
        print(f"  Loaded {len(vcf_variants)} variants from known SNPs VCF")
        if sample_names:
            print(f"  Samples in VCF: {len(sample_names)}")
        
        return vcf_variants
    
    except FileNotFoundError:
        print(f"  Warning: Known SNPs VCF file not found: {vcf_file}")
        return None
    except Exception as e:
        print(f"  Warning: Error loading known SNPs VCF file: {e}")
        import traceback
        traceback.print_exc()
        return None


def check_known_snp_overlap(chrom: str, pos: int, ref: str, alt: str, 
                            vcf_variants: dict) -> dict:
    """
    Check if a variant overlaps with one in the known SNPs VCF.
    Returns a dictionary with annotation fields.
    """
    result = {
        'known_snp_overlap': 'No',
        'known_snp_ref': '',
        'known_snp_alt': '',
        'known_snp_qual': '',
        'known_snp_af': '',
        'known_snp_dp': '',
        'known_snp_type': '',
        'known_snp_ns': '',
        'known_snp_match': ''
    }
    
    if not vcf_variants:
        return result
    
    key = (chrom, pos)
    
    if key in vcf_variants:
        variant = vcf_variants[key]
        result['known_snp_overlap'] = 'Yes'
        result['known_snp_ref'] = variant['ref']
        result['known_snp_alt'] = variant['alt']
        result['known_snp_qual'] = variant['qual']
        result['known_snp_af'] = variant['af']
        result['known_snp_dp'] = variant['dp']
        result['known_snp_type'] = variant['type']
        result['known_snp_ns'] = variant['ns']
        
        # Determine match type
        vcf_alts = variant['alt'].split(',')
        if ref == variant['ref'] and alt in vcf_alts:
            result['known_snp_match'] = 'exact'
        elif ref == variant['ref']:
            result['known_snp_match'] = 'same_position_same_ref'
        else:
            result['known_snp_match'] = 'same_position_diff_ref'
    
    return result


def annotate_known_snps(df: pd.DataFrame, vcf_variants: dict) -> pd.DataFrame:
    """Add known SNP annotation columns to dataframe based on coordinate overlap."""
    if vcf_variants is None:
        return df
    
    print("Annotating known SNP overlaps...")
    
    # Apply the overlap check to each row
    annotations = df.apply(
        lambda row: pd.Series(check_known_snp_overlap(
            str(row['CHROM']),
            int(row['POS']),
            str(row['REF']),
            str(row['ALT']),
            vcf_variants
        )),
        axis=1
    )
    
    # Concatenate the new columns
    df = pd.concat([df, annotations], axis=1)
    
    overlap_count = (df['known_snp_overlap'] == 'Yes').sum()
    exact_match = (df['known_snp_match'] == 'exact').sum()
    print(f"  Found {overlap_count} variants at known SNP positions")
    print(f"  Of which {exact_match} are exact matches (same ref and alt)")
    
    return df


# =============================================================================
# Main Processing
# =============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Process VCF files with optional lncRNA and known SNP annotations"
    )
    ap.add_argument("input_vcf", help="Input .vcf or .vcf.gz")
    ap.add_argument("output_csv", help="Output CSV path")
    ap.add_argument("--lncrna", help="lncRNA reference file (GTF-like format)", default=None)
    ap.add_argument("--known_snps", help="Known SNPs/consensus VCF file", default=None)
    args = ap.parse_args()

    print(f"\n{'='*60}")
    print("VCF Processing Pipeline")
    print(f"{'='*60}")
    print(f"Input VCF: {args.input_vcf}")
    print(f"Output CSV: {args.output_csv}")
    if args.lncrna:
        print(f"lncRNA reference: {args.lncrna}")
    if args.known_snps:
        print(f"Known SNPs VCF: {args.known_snps}")
    print(f"{'='*60}\n")

    # Load reference files first
    lncrna_regions = load_lncrna_reference(args.lncrna) if args.lncrna else None
    known_snps = load_known_snps_vcf(args.known_snps) if args.known_snps else None

    # Read and process input VCF
    print(f"\nReading input VCF: {args.input_vcf}")
    df = read_vcf_to_df(args.input_vcf)
    print(f"  Loaded {len(df)} variants")

    # Split INFO column
    if "INFO" in df.columns:
        df = split_info_column(df)
        
        # Parse SnpEff annotations if present
        if "ANN" in df.columns:
            print("Parsing SnpEff annotations...")
            ann_df = df["ANN"].apply(parse_snpeff_ann_first)
            df = pd.concat([df, ann_df], axis=1)

    # Parse FORMAT/sample columns
    if all(c in df.columns for c in ["FORMAT", "NORMAL", "TUMOR"]):
        print("Parsing sample genotype data...")
        format_df = df.apply(
            lambda row: parse_format_row(row["FORMAT"], row["NORMAL"], row["TUMOR"]),
            axis=1
        )
        df = pd.concat([df.drop(columns=["FORMAT", "NORMAL", "TUMOR"]), format_df], axis=1)

    # Add lncRNA annotations
    if lncrna_regions:
        df = annotate_lncrna(df, lncrna_regions)

    # Add known SNP annotations
    if known_snps:
        df = annotate_known_snps(df, known_snps)

    # Drop unwanted columns
    drop_cols = ["ID", "QUAL", "GPV", "FILTER", "SOMATIC", "GQ_normal", "GQ_tumor"]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")

    # Write output
    print(f"\nWriting output to: {args.output_csv}")
    df.to_csv(args.output_csv, index=False)
    print(f"  Wrote {len(df)} rows with {len(df.columns)} columns")
    
    # Print column summary
    print(f"\nOutput columns: {', '.join(df.columns[:10])}...")
    if 'lncRNA' in df.columns:
        print(f"  lncRNA column included")
    if 'known_snp_overlap' in df.columns:
        print(f"  Known SNP columns included")
    
    print(f"\n{'='*60}")
    print("Processing complete!")
    print(f"{'='*60}\n")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
