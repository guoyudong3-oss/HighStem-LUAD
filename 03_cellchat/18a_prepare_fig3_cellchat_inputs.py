#!/usr/bin/env python3
"""
Prepare Fig3 CellChat inputs from full integrated h5ad + latest malignant stem labels.

Design goals:
- Follow paper-like analysis flow (CellChat centered), not paper result matching.
- Save step-wise outputs so failures can resume without full rerun.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Dict, Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare checkpointed CellChat inputs for Fig3 analysis."
    )
    p.add_argument("--run-dir", required=True, help="Run directory for outputs.")
    p.add_argument(
        "--integrated-h5ad",
        default="data/processed/discovery_integrated.h5ad",
        help="Integrated full-cell h5ad.",
    )
    p.add_argument(
        "--malignant-map-csv",
        required=True,
        help="CSV with columns: cell_id, stem_group (from latest malignant set).",
    )
    p.add_argument(
        "--max-per-group",
        type=int,
        default=1500,
        help="Downsample cap per cellchat_label (0 means no cap).",
    )
    p.add_argument(
        "--min-cells-per-group",
        type=int,
        default=30,
        help="Drop labels with fewer cells than this threshold.",
    )
    p.add_argument(
        "--celltype-annot-col",
        default="celltype_annot",
        help="Preferred column name for original TME annotation label.",
    )
    p.add_argument(
        "--nonmalignant-label-mode",
        choices=["compact", "original"],
        default="compact",
        help="How to label non-malignant cells in CellChat groups.",
    )
    p.add_argument(
        "--sample-id",
        default="",
        help="Optional sample ID for per-sample CellChat prep (empty = all samples).",
    )
    p.add_argument(
        "--sample-col",
        default="sample_id",
        help="Preferred sample column in obs when --sample-id is set.",
    )
    p.add_argument(
        "--expr-source",
        choices=["x", "raw_log1p"],
        default="x",
        help=(
            "Expression source exported to CellChat. "
            "'x' uses adata.X directly; 'raw_log1p' uses adata.raw counts "
            "with per-cell normalization (scale-factor) + log1p."
        ),
    )
    p.add_argument(
        "--raw-gene-mode",
        choices=["matched_var", "all_raw"],
        default="matched_var",
        help=(
            "When --expr-source raw_log1p: 'matched_var' keeps only genes in adata.var_names "
            "(usually current HVGs); 'all_raw' exports all raw genes."
        ),
    )
    p.add_argument(
        "--raw-scale-factor",
        type=float,
        default=10000.0,
        help="Scale factor for raw_log1p normalization (default 1e4).",
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_json(obj: dict, path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)


def normalize_label(x: object) -> str:
    if x is None:
        return ""
    s = str(x).strip()
    if s.lower() in {"nan", "none", ""}:
        return ""
    return s


def derive_celltype_annot(row: pd.Series, preferred_col: str) -> str:
    for col in (preferred_col, "Cell_type.refined", "Cell_type"):
        val = normalize_label(row.get(col, ""))
        if val:
            return val
    return ""


def map_tme_label(
    cell_id: str,
    row: pd.Series,
    mal_map: Dict[str, str],
    nonmalignant_label_mode: str,
) -> str:
    if cell_id in mal_map:
        return mal_map[cell_id]

    annot = normalize_label(row.get("celltype_annot", ""))
    ct = normalize_label(row.get("Cell_type", ""))
    ctr = normalize_label(row.get("Cell_type.refined", ""))

    if annot == "Epithelial cells" or ctr == "Epithelial cells" or ct == "Epithelial cells":
        return "NormalEpi"

    if nonmalignant_label_mode == "original":
        return annot

    # Keep a compact, robust label set for CellChat.
    remap = {
        "T/NK cells": "T_NK",
        "T lymphocytes": "T",
        "NK cells": "NK",
        "B lymphocytes": "B",
        "Myeloid cells": "Myeloid",
        "Fibroblasts": "Fib",
        "Endothelial cells": "Endo",
        "MAST cells": "Mast",
    }
    for key in (annot, ct, ctr):
        if key in remap:
            return remap[key]

    return ""


def downsample_by_group(
    meta: pd.DataFrame,
    label_col: str,
    max_per_group: int,
    rng: np.random.Generator,
) -> np.ndarray:
    picked: list[str] = []
    for g, sub in meta.groupby(label_col):
        ids = sub.index.to_numpy()
        if max_per_group > 0 and len(ids) > max_per_group:
            ids = rng.choice(ids, size=max_per_group, replace=False)
        picked.extend(ids.tolist())
    return np.array(picked, dtype=object)


def write_provenance(path: Path, lines: Iterable[str]) -> None:
    with path.open("w", encoding="utf-8") as f:
        for ln in lines:
            f.write(str(ln).rstrip() + "\n")


def main() -> int:
    args = parse_args()
    run_dir = ensure_dir(Path(args.run_dir).resolve())
    input_dir = ensure_dir(run_dir / "cellchat_input")
    ckpt_dir = ensure_dir(run_dir / "checkpoints_fig3_prep")
    tables_dir = ensure_dir(run_dir / "tables")

    rng = np.random.default_rng(args.seed)

    paths = {
        "meta_full": input_dir / "cellchat_metadata_full.csv",
        "counts_full": input_dir / "cellchat_label_counts_full.csv",
        "meta_sampled": input_dir / "cellchat_metadata_sampled.csv",
        "counts_sampled": input_dir / "cellchat_label_counts_sampled.csv",
        "expr_gene_by_cell": input_dir / "cellchat_expr_gene_by_cell.csv.gz",
        "manifest": input_dir / "cellchat_input_manifest.json",
        "provenance": tables_dir / "fig3_cellchat_prep_provenance.txt",
    }

    print("=== Fig3 input prep ===")
    print("run_dir:", run_dir)
    print("integrated_h5ad:", args.integrated_h5ad)
    print("malignant_map_csv:", args.malignant_map_csv)
    print("expr_source:", args.expr_source)
    if args.expr_source == "raw_log1p":
        print("raw_gene_mode:", args.raw_gene_mode, "raw_scale_factor:", args.raw_scale_factor)
    if args.sample_id:
        print("sample filter:", args.sample_id, f"(col={args.sample_col})")

    # Step 1: load map
    ck1 = ckpt_dir / "01_malignant_map_loaded.ok"
    map_df = pd.read_csv(args.malignant_map_csv)
    if "cell_id" not in map_df.columns or "stem_group" not in map_df.columns:
        raise ValueError("malignant_map_csv must contain columns: cell_id, stem_group")
    map_df = map_df.dropna(subset=["cell_id", "stem_group"]).copy()
    map_df["cell_id"] = map_df["cell_id"].astype(str)
    map_df["stem_group"] = map_df["stem_group"].astype(str)
    mal_map = dict(zip(map_df["cell_id"], map_df["stem_group"]))
    ck1.write_text("ok\n", encoding="utf-8")
    print(f"[01] malignant map loaded: {len(mal_map)} cells")

    # Step 2: load full h5ad
    ck2 = ckpt_dir / "02_integrated_loaded.ok"
    adata = ad.read_h5ad(args.integrated_h5ad)
    ck2.write_text("ok\n", encoding="utf-8")
    print(f"[02] integrated loaded: cells={adata.n_obs}, genes={adata.n_vars}")

    # Step 2b: optional per-sample subset.
    sample_col_used = ""
    if args.sample_id:
        ck2b = ckpt_dir / "02b_sample_subset.ok"
        obs_cols = list(adata.obs.columns)
        for cand in (args.sample_col, "sample_id", "Sample"):
            if cand in obs_cols:
                sample_col_used = cand
                break
        if not sample_col_used:
            raise ValueError(
                "--sample-id is set but no sample column found in obs. "
                f"Tried: {args.sample_col}, sample_id, Sample."
            )
        sid = str(args.sample_id)
        sample_vec = adata.obs[sample_col_used].astype(str).to_numpy()
        keep_mask = sample_vec == sid
        n_keep = int(np.sum(keep_mask))
        if n_keep == 0:
            raise ValueError(
                f"Sample '{sid}' not found in column '{sample_col_used}'."
            )
        adata = adata[keep_mask, :].copy()
        ck2b.write_text("ok\n", encoding="utf-8")
        print(
            f"[02b] sample subset: col={sample_col_used}, sample={sid}, "
            f"kept cells={adata.n_obs}"
        )

    # Step 3: build labels
    ck3 = ckpt_dir / "03_labels_built.ok"
    obs = adata.obs.copy()
    obs["celltype_annot"] = [
        derive_celltype_annot(obs.loc[cid], args.celltype_annot_col)
        for cid in obs.index.astype(str)
    ]
    labels = [
        map_tme_label(
            str(cid),
            obs.loc[cid],
            mal_map,
            args.nonmalignant_label_mode,
        )
        for cid in obs.index.astype(str)
    ]
    obs["cellchat_label"] = labels
    obs["is_new_malignant"] = obs.index.astype(str).isin(mal_map)
    obs = obs[obs["cellchat_label"] != ""].copy()
    counts_full = (
        obs["cellchat_label"]
        .value_counts(dropna=False)
        .rename_axis("cellchat_label")
        .reset_index(name="n_cells")
    )
    if args.min_cells_per_group > 1:
        keep = set(
            counts_full.loc[counts_full["n_cells"] >= args.min_cells_per_group, "cellchat_label"]
            .astype(str)
            .tolist()
        )
        obs = obs[obs["cellchat_label"].astype(str).isin(keep)].copy()
        counts_full = (
            obs["cellchat_label"]
            .value_counts(dropna=False)
            .rename_axis("cellchat_label")
            .reset_index(name="n_cells")
        )
    obs.to_csv(paths["meta_full"])
    counts_full.to_csv(paths["counts_full"], index=False)
    ck3.write_text("ok\n", encoding="utf-8")
    print(f"[03] labels built: kept cells={obs.shape[0]}, groups={obs['cellchat_label'].nunique()}")

    # Step 4: downsample
    ck4 = ckpt_dir / "04_downsampled.ok"
    sampled_ids = downsample_by_group(
        meta=obs,
        label_col="cellchat_label",
        max_per_group=args.max_per_group,
        rng=rng,
    )
    obs_sampled = obs.loc[sampled_ids].copy()
    counts_sampled = (
        obs_sampled["cellchat_label"]
        .value_counts(dropna=False)
        .rename_axis("cellchat_label")
        .reset_index(name="n_cells")
    )
    obs_sampled.to_csv(paths["meta_sampled"])
    counts_sampled.to_csv(paths["counts_sampled"], index=False)
    ck4.write_text("ok\n", encoding="utf-8")
    print(
        f"[04] downsampled: cells={obs_sampled.shape[0]}, groups={obs_sampled['cellchat_label'].nunique()}"
    )

    # Step 5: export expression (gene x cell) for sampled cells.
    ck5 = ckpt_dir / "05_expr_exported.ok"
    ad_sub = adata[obs_sampled.index.astype(str).to_numpy(), :].copy()
    if args.expr_source == "x":
        X = ad_sub.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float32)
        expr_df = pd.DataFrame(
            X,
            index=ad_sub.obs_names.astype(str),
            columns=ad_sub.var_names.astype(str),
        ).T
    else:
        if adata.raw is None:
            raise ValueError(
                "--expr-source raw_log1p requested, but adata.raw is missing."
            )
        raw_all = adata.raw[obs_sampled.index.astype(str).to_numpy(), :]
        if args.raw_gene_mode == "matched_var":
            var_main = pd.Index(ad_sub.var_names.astype(str))
            var_raw = pd.Index(raw_all.var_names.astype(str))
            keep_genes = var_main.intersection(var_raw)
            if len(keep_genes) == 0:
                raise ValueError(
                    "No overlapping genes between adata.var_names and adata.raw.var_names "
                    "for --raw-gene-mode matched_var."
                )
            raw_all = raw_all[:, keep_genes.to_numpy()]
            gene_names = keep_genes.to_numpy().astype(str)
        else:
            gene_names = raw_all.var_names.astype(str)

        X_raw = raw_all.X
        if sparse.issparse(X_raw):
            X_raw = X_raw.tocsr().astype(np.float32, copy=False)
            libsize = np.asarray(X_raw.sum(axis=1)).reshape(-1).astype(np.float32)
            libsize[libsize <= 0] = 1.0
            scale = (float(args.raw_scale_factor) / libsize).astype(np.float32)
            X_norm = X_raw.multiply(scale[:, None])
            X_norm.data = np.log1p(X_norm.data).astype(np.float32, copy=False)
            X = X_norm.toarray()
        else:
            X = np.asarray(X_raw, dtype=np.float32)
            libsize = X.sum(axis=1, dtype=np.float64)
            libsize[libsize <= 0] = 1.0
            X *= (float(args.raw_scale_factor) / libsize[:, None]).astype(np.float32)
            np.log1p(X, out=X)

        expr_df = pd.DataFrame(
            X,
            index=ad_sub.obs_names.astype(str),
            columns=gene_names,
        ).T
    expr_df.to_csv(paths["expr_gene_by_cell"], compression="gzip")
    ck5.write_text("ok\n", encoding="utf-8")
    print(f"[05] expression exported: {expr_df.shape[0]} genes x {expr_df.shape[1]} cells")

    # Step 6: manifest + provenance
    manifest = {
        "run_dir": str(run_dir),
        "integrated_h5ad": args.integrated_h5ad,
        "malignant_map_csv": args.malignant_map_csv,
        "seed": args.seed,
        "max_per_group": args.max_per_group,
        "min_cells_per_group": args.min_cells_per_group,
        "celltype_annot_col": args.celltype_annot_col,
        "nonmalignant_label_mode": args.nonmalignant_label_mode,
        "sample_id": args.sample_id,
        "sample_col": args.sample_col,
        "sample_col_used": sample_col_used,
        "expr_source": args.expr_source,
        "raw_gene_mode": args.raw_gene_mode,
        "raw_scale_factor": args.raw_scale_factor,
        "n_integrated_cells": int(adata.n_obs),
        "n_integrated_genes": int(adata.n_vars),
        "n_malignant_map_cells": int(len(mal_map)),
        "n_labeled_cells_full": int(obs.shape[0]),
        "n_labeled_groups_full": int(obs["cellchat_label"].nunique()),
        "n_labeled_cells_sampled": int(obs_sampled.shape[0]),
        "n_labeled_groups_sampled": int(obs_sampled["cellchat_label"].nunique()),
        "paths": {k: str(v) for k, v in paths.items()},
    }
    save_json(manifest, paths["manifest"])
    write_provenance(
        paths["provenance"],
        [
            "Fig3 CellChat input preparation",
            f"Integrated h5ad: {args.integrated_h5ad}",
            f"Malignant map: {args.malignant_map_csv}",
            f"Mapped malignant labels: {len(mal_map)}",
            f"Final sampled cells: {obs_sampled.shape[0]}",
            f"Final groups: {obs_sampled['cellchat_label'].nunique()}",
            f"Groups: {', '.join(counts_sampled['cellchat_label'].astype(str).tolist())}",
            (
                "Method note: labels are derived from latest malignant stem groups + "
                f"TME annotation mode={args.nonmalignant_label_mode} "
                f"(preferred column={args.celltype_annot_col}); "
                f"expr_source={args.expr_source}, raw_gene_mode={args.raw_gene_mode}, "
                f"raw_scale_factor={args.raw_scale_factor}."
            ),
            (
                "Sample note: "
                + (
                    f"per-sample run for {args.sample_id} "
                    f"(obs column={sample_col_used or args.sample_col})"
                    if args.sample_id
                    else "no sample filter (all samples)."
                )
            ),
        ],
    )
    print("[06] manifest/provenance saved")

    print("✅ Done. Files:")
    for k, v in paths.items():
        print(f" - {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
