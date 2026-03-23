#!/usr/bin/env python3

import argparse
import os
from typing import Dict, List, Optional

import pandas as pd


def infer_source(filepath: str) -> str:
    base = os.path.basename(filepath)
    for suffix in (".processed.csv", ".csv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def format_df(df: pd.DataFrame) -> pd.DataFrame:
    for col in ("FREQ_normal", "FREQ_tumor"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.rstrip("%").astype(float)

    if "DP4_normal" in df.columns:
        df["Total_depth_normal"] = (
            df["DP4_normal"].astype(str).str.split(",")
            .apply(lambda x: sum(map(int, x)))
        )
    if "DP4_tumor" in df.columns:
        df["Total_depth_exp"] = (
            df["DP4_tumor"].astype(str).str.split(",")
            .apply(lambda x: sum(map(int, x)))
        )

    df.rename(columns={"FREQ_tumor": "FREQ_exp"}, inplace=True)

    df.drop(columns=["ANN"], inplace=True, errors="ignore")

    block = [
        "FREQ_normal", "FREQ_exp",
        "DP4_normal", "Total_depth_normal",
        "DP4_tumor", "Total_depth_exp",
        "gene_id", "source", "shared_count", "Total Count",
    ]

    lncrna_cols = ["lncRNA"]
    known_snp_cols = [
        "known_snp_overlap", "known_snp_ref", "known_snp_alt",
        "known_snp_qual", "known_snp_af", "known_snp_dp",
        "known_snp_type", "known_snp_ns", "known_snp_match",
    ]

    cols = list(df.columns)
    block_present = [c for c in block if c in cols]
    lncrna_present = [c for c in lncrna_cols if c in cols]
    known_snp_present = [c for c in known_snp_cols if c in cols]
    special_cols = set(block_present + lncrna_present + known_snp_present)
    cols_wo_special = [c for c in cols if c not in special_cols]

    if "SPV" in cols_wo_special and "LOF" in cols_wo_special:
        spv_i = cols_wo_special.index("SPV")
        cols_new = (
            cols_wo_special[: spv_i + 1]
            + block_present
            + cols_wo_special[spv_i + 1 :]
            + lncrna_present
            + known_snp_present
        )
        df = df[cols_new]
    elif lncrna_present or known_snp_present:
        other_cols = [c for c in cols if c not in special_cols]
        df = df[other_cols + block_present + lncrna_present + known_snp_present]

    return df


def combine_dataframes(dataframes: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    if not dataframes:
        return pd.DataFrame()

    df_list = [
        df.assign(source=source).drop(columns=["shared_count"], errors="ignore")
        for source, df in dataframes.items()
    ]
    combined = pd.concat(df_list, ignore_index=True)

    if not {"CHROM", "POS"}.issubset(combined.columns):
        raise ValueError("Missing required columns CHROM and/or POS")

    combined["variant_id"] = combined["CHROM"].astype(str) + "_" + combined["POS"].astype(str)
    shared_counts = combined.groupby("variant_id")["source"].nunique().rename("shared_count")
    combined = combined.merge(shared_counts, on="variant_id", how="left")
    
    return combined.drop(columns=["variant_id"]).reset_index(drop=True)


def read_group_files(files: Optional[List[str]], index_col: Optional[int] = None) -> Dict[str, pd.DataFrame]:
    if not files:
        return {}
    return {infer_source(fp): pd.read_csv(fp, index_col=index_col) for fp in files}


def write_output(dataframes: Dict[str, pd.DataFrame], out_path: str, drop_cols: List[str]) -> None:
    if not dataframes:
        pd.DataFrame().to_csv(out_path, index=False)
        return

    combined = combine_dataframes(dataframes)
    
    combined.drop(columns=[c for c in drop_cols if c in combined.columns], inplace=True)
    
    combined["Total Count"] = len(dataframes)

    combined = format_df(combined)
    combined.to_csv(out_path, index=False)


def main():
    p = argparse.ArgumentParser(
        description="Combine per-sample processed CSVs into per-group SNP/indel outputs."
    )
    p.add_argument("--group", required=True, help="Group name for output filenames")
    p.add_argument("--snps", nargs="*", default=[], help="SNP processed CSV files")
    p.add_argument("--indels", nargs="*", default=[], help="Indel processed CSV files")
    p.add_argument("-o", "--outdir", default=".", help="Output directory")
    p.add_argument("--index-col", type=int, default=None, help="Index column for read_csv")
    p.add_argument(
        "--drop-cols", nargs="*",
        default=[
            "GT_normal", "GT_tumor", "DP_normal", "DP_tumor",
            "RD_normal", "RD_tumor", "AD_normal", "AD_tumor", "DP"
        ],
        help="Columns to drop if present"
    )
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for variant_type, files in [("snps", args.snps), ("indels", args.indels)]:
        dataframes = read_group_files(files, args.index_col)
        out_path = os.path.join(args.outdir, f"{args.group}_{variant_type}.csv")
        write_output(dataframes, out_path, args.drop_cols)


if __name__ == "__main__":
    main()
