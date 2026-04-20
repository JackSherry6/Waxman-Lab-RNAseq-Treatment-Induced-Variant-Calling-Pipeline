#!/usr/bin/env python3
import re
import gzip
import argparse
import pandas as pd
from collections import defaultdict


SNPEFF_KEYS = [
    "snpeff_allele", "snpeff_effect", "snpeff_impact",
    "snpeff_gene_name", "snpeff_gene_id", "snpeff_feature_type",
    "snpeff_feature_id", "snpeff_biotype", "snpeff_rank",
    "snpeff_hgvs_c", "snpeff_hgvs_p", "snpeff_cdna_poslen",
    "snpeff_cds_poslen", "snpeff_aa_poslen", "snpeff_distance",
    "snpeff_errors_warnings_info",
]


def parse_snpeff_ann_all(ann_value) -> list:
    empty = {k: None for k in SNPEFF_KEYS}

    if pd.isna(ann_value):
        return [empty]

    ann_str = str(ann_value).strip()
    if not ann_str or ann_str == ".":
        return [empty]

    entries = ann_str.split(",")
    parsed = []
    for entry in entries:
        fields = entry.split("|")
        fields += [""] * (len(SNPEFF_KEYS) - len(fields))
        row = {}
        for i, key in enumerate(SNPEFF_KEYS):
            row[key] = fields[i] if fields[i] else None
        parsed.append(row)

    return parsed if parsed else [empty]


def expand_ann_column(df: pd.DataFrame) -> pd.DataFrame:
    if "ANN" not in df.columns:
        return df

    print("Expanding SnpEff annotations (all entries per variant)...")

    rows = []
    for idx, row in df.iterrows():
        ann_entries = parse_snpeff_ann_all(row["ANN"])
        row_dict = row.to_dict()
        for ann in ann_entries:
            new_row = {**row_dict, **ann}
            rows.append(new_row)

    expanded = pd.DataFrame(rows)
    n_orig = len(df)
    n_new = len(expanded)
    print(f"  Expanded {n_orig} rows -> {n_new} rows "
          f"({n_new - n_orig} additional from multi-gene annotations)")
    return expanded

def compute_depth_from_dp4(dp4_value):
    if pd.isna(dp4_value):
        return None
    try:
        return sum(int(x) for x in str(dp4_value).split(","))
    except:
        return None

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
        print(f"  Loaded {total_regions} lncRNA exon regions across "
              f"{len(lncrna_regions)} chromosomes")

        return dict(lncrna_regions)

    except FileNotFoundError:
        print(f"  Warning: lncRNA file not found: {lncrna_file}")
        return None
    except Exception as e:
        print(f"  Warning: Error loading lncRNA file: {e}")
        return None


def check_lncrna_overlap(chrom: str, pos: int, lncrna_regions: dict) -> str:
    
    if not lncrna_regions or chrom not in lncrna_regions:
        return ""

    overlapping_genes = []

    for region in lncrna_regions[chrom]:
        if region['end'] < pos:
            continue
        if region['start'] > pos:
            break
        overlapping_genes.append(region['gene_name'])

    unique_genes = sorted(dict.fromkeys(overlapping_genes))
    return ';'.join(unique_genes) if unique_genes else ""


def annotate_lncrna(df: pd.DataFrame, lncrna_regions: dict) -> pd.DataFrame:
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

def load_known_snps_vcf(vcf_file: str) -> dict:
    
    if not vcf_file:
        return None

    print(f"Loading known SNPs VCF from {vcf_file}...")

    vcf_variants = {}
    sample_names = []

    try:
        opener = gzip.open if vcf_file.endswith('.gz') else open

        with opener(vcf_file, 'rt') as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith('#CHROM'):
                    parts = line.split('\t')
                    if len(parts) > 9:
                        sample_names = parts[9:]
                    continue

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

                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[item] = 'True'

                samples_with_variant = 0
                sample_genotypes = []
                if len(parts) > 9:
                    for sample_data in parts[9:]:
                        gt = sample_data.split(':')[0] if sample_data else '.'
                        if gt not in ['.', './.', '.|.']:
                            sample_genotypes.append(gt)
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
                    'ns': info_dict.get('NS', ''),
                    'ac': info_dict.get('AC', ''),
                    'an': info_dict.get('AN', ''),
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
        return df

    print("Annotating known SNP overlaps...")

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
    print(f"  Found {overlap_count} variants at known SNP positions")
    print(f"  Of which {exact_match} are exact matches (same ref and alt)")

    return df

def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [
        "GT_normal", "GT_tumor",
        "DP_normal", "DP_tumor",
        "RD_normal", "RD_tumor",
        "AD_normal", "AD_tumor",
        "FREQ_normal", "FREQ_tumor",
        "DP4_normal", "DP4_tumor",
    ]

    sample_cols_present = [c for c in sample_cols if c in df.columns]

    if not sample_cols_present:
        print("  No sample columns found to reorder")
        return df

    if "Depth_normal" not in df.columns:
        df["Depth_normal"] = ""
    if "Depth_exp" not in df.columns:
        df["Depth_exp"] = ""

    insert_order = []
    for c in sample_cols_present:
        insert_order.append(c)
        if c == "DP4_normal":
            insert_order.append("Depth_normal")
        elif c == "DP4_tumor":
            insert_order.append("Depth_exp")

    if "Depth_normal" not in insert_order:
        insert_order.append("Depth_normal")
    if "Depth_exp" not in insert_order:
        insert_order.append("Depth_exp")

    cols_to_move = set(insert_order)
    remaining = [c for c in df.columns if c not in cols_to_move]

    if "SPV" in remaining:
        spv_idx = remaining.index("SPV") + 1
        new_order = remaining[:spv_idx] + insert_order + remaining[spv_idx:]
    elif "ANN" in remaining:
        ann_idx = remaining.index("ANN")
        new_order = remaining[:ann_idx] + insert_order + remaining[ann_idx:]
    else:
        new_order = remaining + insert_order

    print("  Reordered sample columns between SPV and ANN")
    return df[new_order]

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

    lncrna_regions = load_lncrna_reference(args.lncrna) if args.lncrna else None
    known_snps = load_known_snps_vcf(args.known_snps) if args.known_snps else None

    print(f"\nReading input VCF: {args.input_vcf}")
    df = read_vcf_to_df(args.input_vcf)
    print(f"  Loaded {len(df)} variants")

    if "INFO" in df.columns:
        df = split_info_column(df)

        if "ANN" in df.columns:
            df = expand_ann_column(df)

    if all(c in df.columns for c in ["FORMAT", "NORMAL", "TUMOR"]):
        print("Parsing sample genotype data...")
        format_df = df.apply(
            lambda row: parse_format_row(row["FORMAT"], row["NORMAL"], row["TUMOR"]),
            axis=1
        )
        df = pd.concat([df.drop(columns=["FORMAT", "NORMAL", "TUMOR"]), format_df], axis=1)

    if "DP4_normal" in df.columns:
        df["Depth_normal"] = df["DP4_normal"].apply(compute_depth_from_dp4)
    if "DP4_tumor" in df.columns:
        df["Depth_exp"] = df["DP4_tumor"].apply(compute_depth_from_dp4)
    if "DP4_exp" in df.columns:
        df["Depth_exp"] = df["DP4_exp"].apply(compute_depth_from_dp4)

    if lncrna_regions:
        df = annotate_lncrna(df, lncrna_regions)

    if known_snps:
        df = annotate_known_snps(df, known_snps)

    print("\nReordering columns...")
    df = reorder_columns(df)

    drop_cols = [
        "ID", "QUAL", "GPV", "FILTER", "SOMATIC", "GQ_normal", "GQ_tumor",
        "DP", "ANN", "snpeff_allele",
        "RD_normal", "AD_normal", "RD_tumor", "AD_tumor","GT_normal","GT_tumor","DP_normal","DP_tumor",
    ]
    print("Dropping unwanted columns...")
    actually_dropped = [c for c in drop_cols if c in df.columns]
    df = df.drop(columns=actually_dropped, errors="ignore")
    if actually_dropped:
        print(f"  Dropped: {', '.join(actually_dropped)}")

    rename_map = {}
    if "DP4_tumor" in df.columns:
        rename_map["DP4_tumor"] = "DP4_exp"
    if "FREQ_tumor" in df.columns:
        rename_map["FREQ_tumor"] = "FREQ_exp"
    if rename_map:
        df = df.rename(columns=rename_map)
        print(f"  Renamed: {rename_map}")

    print(f"\nWriting output to: {args.output_csv}")
    df.to_csv(args.output_csv, index=False)
    print(f"  Wrote {len(df)} rows with {len(df.columns)} columns")

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
