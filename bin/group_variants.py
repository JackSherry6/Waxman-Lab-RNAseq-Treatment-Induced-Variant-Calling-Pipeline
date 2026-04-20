#!/usr/bin/env python3
import argparse
import os
import re
from typing import Dict, List, Optional
import pandas as pd

def infer_source(filepath: str) -> str:
    base = os.path.basename(filepath)

    for suffix in (".processed.csv", ".csv"):
        if base.endswith(suffix):
            stem = base[: -len(suffix)]
            break
    else:
        stem = os.path.splitext(base)[0]

    m = re.search(r"(\d+)_(?:snp|indel)$", stem, re.IGNORECASE)
    if m:
        return m.group(1)

    m = re.search(r"(\d+)$", stem)
    return m.group(1) if m else stem

FREQ_COLS  = ["FREQ_normal", "FREQ_exp"]
DEPTH_COLS = ["Depth_normal", "Depth_exp"]
STAT_COLS  = ["SPV", "SSC"]   # averaged

AVG_COLS = FREQ_COLS + DEPTH_COLS + STAT_COLS


# ── Columns joined across samples ────────────────────────────────────────────
JOIN_COLS = [
    "GT_normal", "GT_tumor",
    "DP_normal", "DP_tumor",
    "DP4_normal", "DP4_exp",
]


SPECIAL_COLS = {"shared_count", "Total Count", "variant_id"}

VARIANT_ID_COLS = [
    "CHROM", "POS", "REF", "ALT",
    "snpeff_effect", "snpeff_impact", "snpeff_gene_name", "snpeff_gene_id",
    "snpeff_feature_type", "snpeff_feature_id", "snpeff_biotype",
    "snpeff_rank", "snpeff_hgvs_c", "snpeff_hgvs_p",
    "snpeff_cdna_poslen", "snpeff_cds_poslen", "snpeff_aa_poslen",
    "snpeff_distance", "snpeff_errors_warnings_info",
    "lncRNA",
    "known_snp_overlap", "known_snp_ref", "known_snp_alt",
    "known_snp_qual", "known_snp_af", "known_snp_dp",
    "known_snp_type", "known_snp_ns", "known_snp_match",
]

def _safe_float(val) -> float:
    try:
        return float(str(val).rstrip("%"))
    except (ValueError, TypeError):
        return float("nan")

def _avg_series(s: pd.Series) -> float:
    vals = s.dropna().apply(_safe_float)
    return vals.mean() if len(vals) else float("nan")

def _join_series(s: pd.Series) -> str:
    return ";".join(s.dropna().astype(str).values)

def aggregate_per_variant(combined: pd.DataFrame) -> pd.DataFrame:
    
    all_cols = list(combined.columns)
    
    group_cols = [c for c in VARIANT_ID_COLS if c in all_cols]
    
    if not group_cols:
        raise ValueError("No variant ID columns found.")

    agg_dict = {}

    avg_present = [c for c in AVG_COLS if c in all_cols]
    join_present = [c for c in JOIN_COLS if c in all_cols]

    for col in avg_present:
        agg_dict[col] = _avg_series

    for col in join_present:
        agg_dict[col] = _join_series


    accounted = set(group_cols) | set(agg_dict.keys()) | SPECIAL_COLS

    for col in all_cols:
        if col not in accounted:
            agg_dict[col] = "first"

    grouped = (
        combined
        .groupby(group_cols, sort=False, dropna=False)
        .agg(agg_dict)
        .reset_index()
    )

    # Rename averages
    rename_dict = {c: f"avg_{c}" for c in avg_present}
    grouped = grouped.rename(columns=rename_dict)


    return grouped


def format_df(df: pd.DataFrame, total_samples: int) -> pd.DataFrame:
    df["total_sample_count"] = total_samples

    if "Total Count" in df.columns:
        df = df.drop(columns=["Total Count"])

    # Block immediately after ALT
    ordered_block = [
        "SS",
        "avg_SPV",
        "avg_SSC",
        "avg_FREQ_normal",
        "avg_FREQ_exp",
        "avg_DP4_normal",
        "avg_Depth_normal",
        "avg_DP4_exp",
        "avg_Depth_exp",
        "shared_count",
        "total_sample_count",
    ]

    ordered_block = [c for c in ordered_block if c in df.columns]

    cols = list(df.columns)

    if "ALT" in cols:
        alt_index = cols.index("ALT")

        remaining = [c for c in cols if c not in ordered_block]

        cols = (
            remaining[:alt_index + 1]
            + ordered_block
            + [c for c in remaining[alt_index + 1:] if c not in ordered_block]
        )

        df = df[cols]

    if "total_sample_count" in df.columns:

        insert_idx = df.columns.get_loc("total_sample_count") + 1

        lof_nmd = [c for c in ["LOF", "NMD"] if c in df.columns]

        for col in lof_nmd:
            df.insert(insert_idx, col, df.pop(col))
            insert_idx += 1

    known_snp_cols = [c for c in df.columns if c.startswith("known_snp")]

    if known_snp_cols:
        remaining = [c for c in df.columns if c not in known_snp_cols]
        df = df[remaining + known_snp_cols]

    return df


def combine_dataframes(dataframes: Dict[str, pd.DataFrame]) -> pd.DataFrame:

    if not dataframes:
        return pd.DataFrame()

    df_list = []

    for source, df in dataframes.items():

        tmp = df.copy()

        for col in ("shared_count", "Total Count", "variant_id"):
            if col in tmp.columns:
                tmp.drop(columns=[col], inplace=True)

        stale = [c for c in tmp.columns if c.startswith("avg_")]
        if stale:
            tmp.drop(columns=stale, inplace=True)

        tmp["source"] = source
        df_list.append(tmp)

    combined = pd.concat(df_list, ignore_index=True)

    if not {"CHROM", "POS"}.issubset(combined.columns):
        raise ValueError("Missing CHROM/POS columns")

    return aggregate_per_variant(combined).reset_index(drop=True)


def read_group_files(
    files: Optional[List[str]], index_col: Optional[int] = None
) -> Dict[str, pd.DataFrame]:

    if not files:
        return {}

    dfs = {}

    for fp in files:

        name = infer_source(fp)

        if name in dfs:
            raise ValueError(f"Duplicate source name: {name}")

        dfs[name] = pd.read_csv(fp, index_col=index_col)

    return dfs


def write_output(
    dataframes: Dict[str, pd.DataFrame],
    out_path: str,
    drop_cols: List[str],
):

    if not dataframes:
        pd.DataFrame().to_csv(out_path, index=False)
        return

    combined = combine_dataframes(dataframes)

    combined.drop(
        columns=[c for c in drop_cols if c in combined.columns],
        inplace=True,
    )

    combined = format_df(combined, len(dataframes))

    combined.to_csv(out_path, index=False)

    print(f"Wrote {len(combined)} rows to {out_path}")


def main():

    p = argparse.ArgumentParser(
        description="Combine per-sample processed CSVs into per-group SNP/indel outputs."
    )

    p.add_argument("--group", required=True)
    p.add_argument("--snps", nargs="*", default=[])
    p.add_argument("--indels", nargs="*", default=[])
    p.add_argument("-o", "--outdir", default=".")
    p.add_argument("--index-col", type=int, default=None)
    p.add_argument("--drop-cols", nargs="*", default=[])

    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for variant_type, files in [
        ("snps", args.snps),
        ("indels", args.indels),
    ]:

        if not files:
            continue

        print(f"\nProcessing {variant_type} ({len(files)} files)...")

        dataframes = read_group_files(files, args.index_col)

        out_path = os.path.join(
            args.outdir,
            f"{args.group}_{variant_type}.csv",
        )

        write_output(dataframes, out_path, args.drop_cols)

    print("\nDone.")


if __name__ == "__main__":
    main()
