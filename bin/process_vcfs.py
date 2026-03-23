#!/usr/bin/env python3
import re
import sys
import gzip
import argparse
import os
import pandas as pd
from collections import defaultdict


def parse_snpeff_ann_first(ann_value: str) -> pd.Series:
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
    keys = str(format_str).split(":")
    row = {}
    for k, n, t in zip(keys, str(normal_str).split(":"), str(tumor_str).split(":")):
        row[f"{k}_normal"] = n
        row[f"{k}_tumor"] = t
    return pd.Series(row)


def parse_gtf_attributes(attr_string: str) -> dict:
    attrs = {}
    pattern = r'(\w+)\s*["\']([^"\']*)["\']'
    matches = re.findall(pattern, attr_string)
    for key, value in matches:
        attrs[key] = value
    return attrs


def load_lncrna_reference(lncrna_file: str) -> dict:
    if not lncrna_file:
        return None
    
    print(f"Loading lncRNA reference from {lncrna_file}...", file=sys.stderr)
    
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
                start = int(parts[3])
                end = int(parts[4])
                attributes = parts[8]
                
                attrs = parse_gtf_attributes(attributes)
                gene_name = attrs.get('gene_name', attrs.get('gene_id', 'unknown'))
                
                lncrna_regions[chrom].append({
                    'start': start,
                    'end': end,
                    'gene_name': gene_name
                })
        
        for chrom in lncrna_regions:
            lncrna_regions[chrom].sort(key=lambda x: x['start'])
        
        total_regions = sum(len(v) for v in lncrna_regions.values())
        print(f"  Loaded {total_regions} lncRNA exon regions across {len(lncrna_regions)} chromosomes", file=sys.stderr)
        
        return dict(lncrna_regions)
    
    except FileNotFoundError:
        print(f"  Warning: lncRNA file not found: {lncrna_file}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"  Warning: Error loading lncRNA file: {e}", file=sys.stderr)
        return None


def check_lncrna_overlap(chrom: str, pos: int, lncrna_regions: dict) -> str:
    if not lncrna_regions or chrom not in lncrna_regions:
        return ""
    
    overlapping_genes = []
    for region in lncrna_regions[chrom]:
        if region['start'] <= pos <= region['end']:
            overlapping_genes.append(region['gene_name'])
        elif region['start'] > pos:
            break
    
    unique_genes = list(dict.fromkeys(overlapping_genes))
    return ';'.join(unique_genes) if unique_genes else ""


def annotate_lncrna(df: pd.DataFrame, lncrna_regions: dict) -> pd.DataFrame:
    if lncrna_regions is None:
        return df
    
    print("Annotating lncRNA overlaps...", file=sys.stderr)
    
    df['lncRNA'] = df.apply(
        lambda row: check_lncrna_overlap(
            str(row['CHROM']),
            int(row['POS']),
            lncrna_regions
        ),
        axis=1
    )
    
    overlap_count = (df['lncRNA'] != '').sum()
    print(f"  Found {overlap_count} variants overlapping lncRNA exons", file=sys.stderr)
    
    return df


def load_known_snps_vcf(vcf_file: str) -> dict:
    if not vcf_file:
        return None
    
    print(f"Loading known SNPs VCF from {vcf_file}...", file=sys.stderr)
    
    vcf_variants = {}
    sample_names = []
    line_count = 0
    
    try:
        if vcf_file.endswith('.gz'):
            f = gzip.open(vcf_file, 'rt')
        else:
            f = open(vcf_file, 'r')
        
        with f:
            for line in f:
                line = line.strip()
                
                if not line:
                    continue
                
                if line.startswith('#CHROM'):
                    parts = line.split('\t')
                    if len(parts) > 9:
                        sample_names = parts[9:]
                    print(f"  Found {len(sample_names)} samples in VCF header", file=sys.stderr)
                    continue
                
                if line.startswith('#'):
                    continue
                
                line_count += 1
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
                
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[item] = 'True'
                
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
                    'ns': info_dict.get('NS', ''),
                    'ac': info_dict.get('AC', ''),
                    'an': info_dict.get('AN', ''),
                }
        
        print(f"  Parsed {line_count} variant lines", file=sys.stderr)
        print(f"  Loaded {len(vcf_variants)} unique positions from known SNPs VCF", file=sys.stderr)
        
        if vcf_variants:
            example_keys = list(vcf_variants.keys())[:3]
            print(f"  Example positions: {example_keys}", file=sys.stderr)
        
        return vcf_variants
    
    except FileNotFoundError:
        print(f"  ERROR: Known SNPs VCF file not found: {vcf_file}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"  ERROR: Failed to load known SNPs VCF: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return None


def check_known_snp_overlap(chrom: str, pos: int, ref: str, alt: str, 
                            vcf_variants: dict) -> dict:
    result = {'known_snp_overlap': 'No',
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
        
        vcf_alts = variant['alt'].split(',')
        if ref == variant['ref'] and alt in vcf_alts:
            result['known_snp_match'] = 'exact'
        elif ref == variant['ref']:
            result['known_snp_match'] = 'same_position_same_ref'
        else:
            result['known_snp_match'] = 'same_position_diff_ref'
    
    return result


def annotate_known_snps(df: pd.DataFrame, vcf_variants: dict) -> pd.DataFrame:
    if vcf_variants is None:
        print("  Skipping known SNP annotation (no variants loaded)", file=sys.stderr)
        return df
    
    print("Annotating known SNP overlaps...", file=sys.stderr)
    
    if len(df) > 0:
        example_positions = list(zip(df['CHROM'].head(3), df['POS'].head(3)))
        print(f"  Example input positions: {example_positions}", file=sys.stderr)
    
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
    
    df = pd.concat([df, annotations], axis=1)
    
    overlap_count = (df['known_snp_overlap'] == 'Yes').sum()
    exact_match = (df['known_snp_match'] == 'exact').sum()
    print(f"  Found {overlap_count} variants at known SNP positions", file=sys.stderr)
    print(f"  Of which {exact_match} are exact matches (same ref and alt)", file=sys.stderr)
    
    return df


def main():
    ap = argparse.ArgumentParser(
        description="Process VCF files with optional lncRNA and known SNP annotations"
    )
    ap.add_argument("input_vcf", help="Input .vcf or .vcf.gz")
    ap.add_argument("output_csv", help="Output CSV path")
    ap.add_argument("--lncrna", help="lncRNA reference file (GTF-like format)", default=None)
    ap.add_argument("--known_snps", help="Known SNPs/consensus VCF file", default=None)
    args = ap.parse_args()

    print(f"\n{'='*60}", file=sys.stderr)
    print("VCF Processing Pipeline", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Input VCF: {args.input_vcf}", file=sys.stderr)
    print(f"Output CSV: {args.output_csv}", file=sys.stderr)
    print(f"lncRNA reference: {args.lncrna}", file=sys.stderr)
    print(f"Known SNPs VCF: {args.known_snps}", file=sys.stderr)
    print(f"{'='*60}\n", file=sys.stderr)

    known_snps_path = args.known_snps
    if known_snps_path:
        basename = os.path.basename(known_snps_path)
        if basename.startswith('NO_FILE') or not os.path.exists(known_snps_path):
            print(f"  Known SNPs file is placeholder or doesn't exist, skipping", file=sys.stderr)
            known_snps_path = None
        elif os.path.getsize(known_snps_path) == 0:
            print(f"  Known SNPs file is empty, skipping", file=sys.stderr)
            known_snps_path = None

    lncrna_path = args.lncrna
    if lncrna_path:
        basename = os.path.basename(lncrna_path)
        if basename.startswith('NO_FILE') or not os.path.exists(lncrna_path):
            print(f"  lncRNA file is placeholder or doesn't exist, skipping", file=sys.stderr)
            lncrna_path = None
        elif os.path.getsize(lncrna_path) == 0:
            print(f"  lncRNA file is empty, skipping", file=sys.stderr)
            lncrna_path = None

    lncrna_regions = load_lncrna_reference(lncrna_path) if lncrna_path else None
    known_snps = load_known_snps_vcf(known_snps_path) if known_snps_path else None

    print(f"Reading input VCF: {args.input_vcf}", file=sys.stderr)
    df = read_vcf_to_df(args.input_vcf)
    print(f"  Loaded {len(df)} variants", file=sys.stderr)

    if "INFO" in df.columns:
        df = split_info_column(df)
        
        if "ANN" in df.columns:
            print("Parsing SnpEff annotations...", file=sys.stderr)
            ann_df = df["ANN"].apply(parse_snpeff_ann_first)
            df = pd.concat([df, ann_df], axis=1)

    if all(c in df.columns for c in ["FORMAT", "NORMAL", "TUMOR"]):
        print("Parsing sample genotype data...", file=sys.stderr)
        format_df = df.apply(
            lambda row: parse_format_row(row["FORMAT"], row["NORMAL"], row["TUMOR"]),
            axis=1
        )
        df = pd.concat([df.drop(columns=["FORMAT", "NORMAL", "TUMOR"]), format_df], axis=1)

    if lncrna_regions:
        df = annotate_lncrna(df, lncrna_regions)

    if known_snps:
        df = annotate_known_snps(df, known_snps)
    else:
        print("  No known SNPs data available for annotation", file=sys.stderr)

    drop_cols = ["ID", "QUAL", "GPV", "FILTER", "SOMATIC", "GQ_normal", "GQ_tumor"]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")

    print(f"\nWriting output to: {args.output_csv}", file=sys.stderr)
    df.to_csv(args.output_csv, index=False)
    print(f"  Wrote {len(df)} rows with {len(df.columns)} columns", file=sys.stderr)
    print(f"  Columns: {list(df.columns)}", file=sys.stderr)
    
    print(f"\n{'='*60}", file=sys.stderr)
    print("Processing complete!", file=sys.stderr)
    print(f"{'='*60}\n", file=sys.stderr)
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())