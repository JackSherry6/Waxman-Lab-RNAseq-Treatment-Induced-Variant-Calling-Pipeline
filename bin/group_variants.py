#!/usr/bin/env python3
"""
Combine per-sample processed CSVs into per-group SNP/indel outputs.
Handles lncRNA and known SNP annotation columns from process_vcfs.py.
"""

import argparse
import os
from typing import Dict, List, Optional

import pandas as pd


def infer_source(filepath: str) -> str:
    """Infer sample/source label from processed CSV filename."""
    base = os.path.basename(filepath)
    for suffix in (".processed.csv", ".csv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def format_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Format and reorder dataframe columns:
    - Convert FREQ columns from percent strings to float
    - Compute total depths from DP4 strings
    - Rename tumor -> exp columns
    - Reorder columns into logical groups
    """
    # Convert percent strings to float
    for col in ("FREQ_normal", "FREQ_tumor"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.rstrip("%").astype(float)

    # Compute total depths from DP4 strings (e.g., "137,64,0,0")
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

    # Rename tumor -> exp for clarity
    df.rename(columns={"FREQ_tumor": "FREQ_exp"}, inplace=True)

    # Drop raw annotation column if present
    df.drop(columns=["ANN"], inplace=True, errors="ignore")

    # Define column ordering blocks
    block = [
        "FREQ_normal", "FREQ_exp",
        "DP4_normal", "Total_depth_normal",
        "DP4_tumor", "Total_depth_exp",
        "gene_id", "source", "shared_count", "Total Count",
    ]

    # lncRNA and known SNP columns to keep together at the end
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
    
    # Remove all special columns from main list
    special_cols = set(block_present + lncrna_present + known_snp_present)
    cols_wo_special = [c for c in cols if c not in special_cols]

    # Insert block between SPV and LOF, put annotation cols at end
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
        # If SPV/LOF not present, just ensure annotation cols are at end
        other_cols = [c for c in cols if c not in special_cols]
        df = df[other_cols + block_present + lncrna_present + known_snp_present]

    return df


def combine_dataframes(dataframes: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Stack sample dataframes and compute shared_count 
    (number of unique samples containing each CHROM+POS).
    """
    if not dataframes:
        return pd.DataFrame()

    # Add source column and stack
    df_list = [
        df.assign(source=source).drop(columns=["shared_count"], errors="ignore")
        for source, df in dataframes.items()
    ]
    combined = pd.concat(df_list, ignore_index=True)

    # Validate required columns
    if not {"CHROM", "POS"}.issubset(combined.columns):
        raise ValueError("Missing required columns CHROM and/or POS")

    # Compute shared_count via groupby
    combined["variant_id"] = combined["CHROM"].astype(str) + "_" + combined["POS"].astype(str)
    shared_counts = combined.groupby("variant_id")["source"].nunique().rename("shared_count")
    combined = combined.merge(shared_counts, on="variant_id", how="left")
    
    return combined.drop(columns=["variant_id"]).reset_index(drop=True)


def read_group_files(files: Optional[List[str]], index_col: Optional[int] = None) -> Dict[str, pd.DataFrame]:
    """Read CSV files into a dictionary keyed by inferred source name."""
    if not files:
        return {}
    return {infer_source(fp): pd.read_csv(fp, index_col=index_col) for fp in files}


def write_output(dataframes: Dict[str, pd.DataFrame], out_path: str, drop_cols: List[str]) -> None:
    """Combine dataframes, format, and write to CSV."""
    if not dataframes:
        pd.DataFrame().to_csv(out_path, index=False)
        return

    combined = combine_dataframes(dataframes)
    
    # Drop unwanted columns
    combined.drop(columns=[c for c in drop_cols if c in combined.columns], inplace=True)
    
    # Add total sample count
    combined["Total Count"] = len(dataframes)

    # Format and reorder
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

    # Process SNPs and indels
    for variant_type, files in [("snps", args.snps), ("indels", args.indels)]:
        dataframes = read_group_files(files, args.index_col)
        out_path = os.path.join(args.outdir, f"{args.group}_{variant_type}.csv")
        write_output(dataframes, out_path, args.drop_cols)


if __name__ == "__main__":
    main()


"""
#!/usr/bin/env python3
import argparse
import os
from typing import Dict, List, Optional

import pandas as pd


def infer_source(filepath: str) -> str:
    base = os.path.basename(filepath)
    for suffix in [".processed.csv", ".csv"]:
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def format_df(df: pd.DataFrame) -> pd.DataFrame:
    # Convert percent strings to float (strip trailing % if present)
    if "FREQ_normal" in df.columns:
        df["FREQ_normal"] = df["FREQ_normal"].astype(str).str.rstrip("%").astype(float)

    if "FREQ_tumor" in df.columns:
        df["FREQ_tumor"] = df["FREQ_tumor"].astype(str).str.rstrip("%").astype(float)

    # Drop ANN if present
    df.drop(columns=["ANN"], inplace=True, errors="ignore")

    # Compute total depths from DP4 strings like "137,64,0,0"
    if "DP4_normal" in df.columns:
        df["Total_depth_normal"] = df["DP4_normal"].astype(str).str.split(",").apply(
            lambda x: sum(map(int, x))
        )

    if "DP4_tumor" in df.columns:
        df["Total_depth_exp"] = df["DP4_tumor"].astype(str).str.split(",").apply(
            lambda x: sum(map(int, x))
        )

    # Rename FREQ_tumor -> FREQ_exp
    if "FREQ_tumor" in df.columns:
        df.rename(columns={"FREQ_tumor": "FREQ_exp"}, inplace=True)

    # Reorder: put this block (in this order) between SPV and LOF
    block = [
        "FREQ_normal", "FREQ_exp",
        "DP4_normal", "Total_depth_normal",
        "DP4_tumor", "Total_depth_exp",
        "gene_id", "source", "shared_count", "Total Count"
    ]

    cols = list(df.columns)
    block_present = [c for c in block if c in cols]
    cols_wo_block = [c for c in cols if c not in block_present]

    if "SPV" in cols_wo_block and "LOF" in cols_wo_block:
        spv_i = cols_wo_block.index("SPV")
        cols_new = cols_wo_block[:spv_i + 1] + block_present + cols_wo_block[spv_i + 1:]
        df = df[cols_new]

    return df


def combine_dataframes_keep_samples(dataframes_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    df_list = []
    for source, df in dataframes_dict.items():
        df_copy = df.copy()
        if "shared_count" in df_copy.columns:
            df_copy = df_copy.drop(columns=["shared_count"])
        df_copy["source"] = source
        df_list.append(df_copy)

    if not df_list:
        return pd.DataFrame()

    combined_df = pd.concat(df_list, ignore_index=True)

    required = {"CHROM", "POS"}
    missing = required - set(combined_df.columns)
    if missing:
        raise ValueError(f"Missing required columns {missing} needed to compute shared_count")

    combined_df["variant_id"] = combined_df["CHROM"].astype(str) + "_" + combined_df["POS"].astype(str)

    shared_counts = (
        combined_df.groupby("variant_id")["source"]
        .nunique()
        .rename("shared_count")
    )

    combined_df = combined_df.merge(shared_counts, on="variant_id", how="left")
    combined_df = combined_df.drop(columns=["variant_id"])
    return combined_df.reset_index(drop=True)


def drop_optional_columns(df: pd.DataFrame, cols_to_drop: List[str]) -> pd.DataFrame:
    existing = [c for c in cols_to_drop if c in df.columns]
    return df.drop(columns=existing) if existing else df


def read_group_files(files: Optional[List[str]], index_col: Optional[int]) -> Dict[str, pd.DataFrame]:
    out: Dict[str, pd.DataFrame] = {}
    if not files:
        return out
    for fp in files:
        source = infer_source(fp)
        out[source] = pd.read_csv(fp, index_col=index_col)
    return out


def write_output(dct: Dict[str, pd.DataFrame], out_path: str, drop_cols: List[str]) -> None:
    if not dct:
        pd.DataFrame().to_csv(out_path, index=False)
        return

    combined = combine_dataframes_keep_samples(dct)
    combined = drop_optional_columns(combined, drop_cols)
    combined["Total Count"] = len(dct)

    # Apply your formatting + reordering after shared_count/Total Count exist
    combined = format_df(combined)

    combined.to_csv(out_path, index=False)


def main():
    p = argparse.ArgumentParser(description="Combine per-sample processed CSVs into per-group SNP/indel outputs.")
    p.add_argument("--group", required=True, help="Group name (used in output filenames).")
    p.add_argument("--snps", nargs="*", default=[], help="SNP processed CSV files for this group.")
    p.add_argument("--indels", nargs="*", default=[], help="Indel processed CSV files for this group.")
    p.add_argument("-o", "--outdir", default=".", help="Output directory (default: current directory).")
    p.add_argument("--index-col", type=int, default=None, help="Index column for pandas.read_csv (default: None).")
    p.add_argument(
        "--drop-cols",
        nargs="*",
        default=["GT_normal", "GT_tumor", "DP_normal", "DP_tumor", "RD_normal", "RD_tumor",
                 "AD_normal", "AD_tumor", "DP"],
        help="Columns to drop if present."
    )
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    snp_dct = read_group_files(args.snps, args.index_col)
    indel_dct = read_group_files(args.indels, args.index_col)

    snp_out = os.path.join(args.outdir, f"{args.group}_snps.csv")
    indel_out = os.path.join(args.outdir, f"{args.group}_indels.csv")

    write_output(snp_dct, snp_out, args.drop_cols)
    write_output(indel_dct, indel_out, args.drop_cols)


if __name__ == "__main__":
    main()

"""