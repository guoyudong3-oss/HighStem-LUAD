#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import zscore


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build Fig5 A-set TCGA activity score tables from run-local high_stem_activity_genes_final.csv."
    )
    p.add_argument("--run-dir", required=True, help="Run directory containing tables/high_stem_activity_genes_final.csv.")
    p.add_argument("--gene-file", default="", help="Optional explicit activity gene CSV.")
    p.add_argument("--expr-file", default="data/external/tcga_luad/expression.tsv.gz")
    p.add_argument("--surv-file", default="data/external/tcga_luad/survival.tsv.gz")
    p.add_argument("--map-file", default="data/external/ensembl_to_symbol.csv")
    p.add_argument("--out-all-csv", default="", help="Output CSV for all samples.")
    p.add_argument("--out-tumor-csv", default="", help="Output CSV for tumor-only samples.")
    return p.parse_args()


def infer_sample_type(sample_id: str) -> str:
    parts = str(sample_id).split("-")
    code = ""
    if len(parts) >= 4:
        code = parts[3][:2]
    if not code and len(sample_id) >= 15:
        code = sample_id[13:15]
    if code == "01":
        return "Tumor"
    if code == "11":
        return "Normal"
    return "Other"


def load_expression(expr_file: Path, map_file: Path) -> pd.DataFrame:
    mapping = pd.read_csv(map_file)
    mapping["gene_id"] = mapping["gene_id"].astype(str).str.split(".").str[0]
    ens_to_sym = pd.Series(mapping.gene_name.values, index=mapping.gene_id).to_dict()

    expr = pd.read_csv(expr_file, sep="\t", index_col=0)
    if str(expr.index[0]).startswith("ENSG"):
        expr.index = expr.index.astype(str).str.split(".").str[0]
        expr.index = expr.index.map(ens_to_sym)
        expr = expr[expr.index.notna()]
        expr = expr.groupby(expr.index).mean()
    return expr


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).resolve()
    gene_file = Path(args.gene_file).resolve() if args.gene_file else run_dir / "tables" / "high_stem_activity_genes_final.csv"
    expr_file = Path(args.expr_file).resolve()
    surv_file = Path(args.surv_file).resolve()
    map_file = Path(args.map_file).resolve()
    out_all_csv = Path(args.out_all_csv).resolve() if args.out_all_csv else run_dir / "tables" / "tcga_highstem_activity_scores_all_fromA.csv"
    out_tumor_csv = Path(args.out_tumor_csv).resolve() if args.out_tumor_csv else run_dir / "tables" / "tcga_highstem_activity_scores_tumor01_fromA.csv"

    if not gene_file.exists():
        raise FileNotFoundError(f"Missing activity gene file: {gene_file}")
    if not expr_file.exists():
        raise FileNotFoundError(f"Missing expression file: {expr_file}")
    if not surv_file.exists():
        raise FileNotFoundError(f"Missing survival file: {surv_file}")
    if not map_file.exists():
        raise FileNotFoundError(f"Missing mapping file: {map_file}")

    genes_df = pd.read_csv(gene_file)
    if "gene" not in genes_df.columns:
        raise ValueError(f"Gene file missing 'gene' column: {gene_file}")
    genes = [str(g) for g in genes_df["gene"].dropna().tolist()]
    if not genes:
        raise ValueError(f"No genes found in {gene_file}")

    expr = load_expression(expr_file, map_file)
    genes_avail = [g for g in genes if g in expr.index]
    if not genes_avail:
        raise ValueError("None of the activity genes were found in the TCGA expression matrix.")

    expr_sub = expr.loc[genes_avail].T.copy()
    expr_sub[genes_avail] = np.log2(expr_sub[genes_avail].clip(lower=0.0) + 1.0)
    expr_z = expr_sub[genes_avail].apply(zscore, axis=0, nan_policy="omit")
    score = expr_z.mean(axis=1)

    score_df = score.to_frame(name="PlasticStem_Score")
    score_df["sample"] = score_df.index.astype(str)
    score_df["sample_short"] = score_df["sample"].str[:15]
    score_df["sample_type"] = score_df["sample"].apply(infer_sample_type)

    surv = pd.read_csv(surv_file, sep="\t")
    surv["sample_short"] = surv["sample"].astype(str).str[:15]
    merged = pd.merge(score_df, surv[["sample_short", "OS", "OS.time"]], on="sample_short", how="left")
    merged = merged[["sample_short", "sample", "PlasticStem_Score", "sample_type", "OS", "OS.time"]].copy()
    tumor_df = merged[merged["sample_type"] == "Tumor"].copy()

    out_all_csv.parent.mkdir(parents=True, exist_ok=True)
    out_tumor_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_all_csv, index=False)
    tumor_df.to_csv(out_tumor_csv, index=False)

    print(f"OUT_ALL={out_all_csv}")
    print(f"OUT_TUMOR={out_tumor_csv}")
    print(f"N_ALL={len(merged)}")
    print(f"N_TUMOR={len(tumor_df)}")
    print(f"N_GENES_USED={len(genes_avail)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
