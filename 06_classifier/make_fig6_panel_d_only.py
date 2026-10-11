#!/usr/bin/env python3
"""
Render only Figure 6 panel D (machine-learning benchmark) with repeated CV.

This is intended for style/method debugging without regenerating the full Figure 6.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List

import anndata as ad
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["pdf.use14corefonts"] = False
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, train_test_split

from make_fig6_reference_like import (
    DEFAULT_RUN_DIR,
    figure6_model_order,
    load_genes,
    model_zoo,
    read_metadata_csv,
    to_dense,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render Figure 6 panel D only.")
    p.add_argument(
        "--run-dir",
        default=DEFAULT_RUN_DIR,
        help="Run directory containing final_consensus_genes.csv and metadata_strict_latest13k.csv",
    )
    p.add_argument(
        "--adata-path",
        default="data/processed/discovery_malignant.h5ad",
        help="Path to atlas-scale h5ad used to recover expression matrix.",
    )
    p.add_argument(
        "--gene-file",
        default="",
        help="Optional explicit gene file (CSV with column 'gene').",
    )
    p.add_argument(
        "--meta-path-classifier",
        default="",
        help="Optional classifier metadata CSV. Defaults to <run-dir>/tables/metadata_strict_latest13k.csv.",
    )
    p.add_argument(
        "--out-dir",
        default="results/reproduce_paper/figures/fig6_panel_d_only",
        help="Directory for panel D outputs.",
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--cv-splits", type=int, default=5)
    p.add_argument("--cv-repeats", type=int, default=10)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def load_classifier_matrix(
    adata_path: str,
    classifier_meta_path: str,
    gene_file: str,
    run_dir: str,
) -> tuple[np.ndarray, np.ndarray, List[str]]:
    adata = ad.read_h5ad(adata_path)
    meta_classifier = read_metadata_csv(classifier_meta_path)
    adata_cells = np.asarray(adata.obs_names.astype(str))
    cls_common = [c for c in adata_cells if c in set(meta_classifier.index.astype(str))]
    if not cls_common:
        raise RuntimeError("No common cells between adata and classifier metadata CSV.")

    adata_classifier = adata[np.asarray(cls_common), :].copy()
    meta_classifier = meta_classifier.loc[cls_common].copy()

    gene_path = gene_file or os.path.join(run_dir, "tables", "final_consensus_genes.csv")
    genes = load_genes(gene_path, adata_classifier)
    if adata_classifier.raw is not None and all(
        g in set(map(str, adata_classifier.raw.var_names)) for g in genes
    ):
        expr_classifier = to_dense(adata_classifier.raw[:, genes].X)
    else:
        expr_classifier = to_dense(adata_classifier[:, genes].X)

    cls_mask = meta_classifier["stem_group_strict"].isin(["HighStem", "LowStem"]).values
    if int(np.sum(cls_mask)) < 50:
        raise RuntimeError("Too few HighStem/LowStem malignant cells for model panels.")

    X = expr_classifier[cls_mask]
    y = (meta_classifier.loc[cls_mask, "stem_group_strict"].values == "HighStem").astype(int)
    return X, y, genes


def repeated_benchmark_df(
    X_train: np.ndarray,
    y_train: np.ndarray,
    seed: int,
    cv_splits: int,
    cv_repeats: int,
) -> pd.DataFrame:
    models = model_zoo(seed)
    rcv = RepeatedStratifiedKFold(
        n_splits=cv_splits,
        n_repeats=cv_repeats,
        random_state=seed,
    )
    recs = []
    for model_name, model in models.items():
        for split_idx, (tr, va) in enumerate(rcv.split(X_train, y_train), start=1):
            fit_model = clone(model)
            fit_model.fit(X_train[tr], y_train[tr])
            prob = fit_model.predict_proba(X_train[va])[:, 1]
            auc = float(roc_auc_score(y_train[va], prob))
            recs.append(
                {
                    "model": model_name,
                    "split_idx": split_idx,
                    "auc": auc,
                }
            )
    return pd.DataFrame(recs)


def d_fill_palette() -> Dict[str, str]:
    return {
        "classif.kknn": "#BEE4D2",
        "classif.lda": "#F4E892",
        "classif.log_reg": "#B9C8FF",
        "classif.naive_bayes": "#F4A2A0",
        "classif.ranger": "#A7D7DB",
        "classif.rpart": "#F6C06C",
        "classif.svm": "#D8F27C",
        "classif.xgboost": "#F5C1CC",
    }


def d_whisker_palette() -> Dict[str, str]:
    return {
        "classif.kknn": "#F26B5B",
        "classif.lda": "#D7B017",
        "classif.log_reg": "#94CF2D",
        "classif.naive_bayes": "#52D88B",
        "classif.ranger": "#11C7CF",
        "classif.rpart": "#2FA9FF",
        "classif.svm": "#C56BFF",
        "classif.xgboost": "#FF4DBA",
    }


def plot_panel_d(
    fit_df: pd.DataFrame,
    out_png: str,
    out_pdf: str,
    dpi: int,
) -> None:
    order = figure6_model_order(sorted(fit_df["model"].unique()))
    grouped = [fit_df.loc[fit_df["model"] == name, "auc"].to_numpy() for name in order]
    fill = d_fill_palette()
    whisk = d_whisker_palette()
    edge = "#4C5BDB"

    fig, ax = plt.subplots(figsize=(13.8, 8.6), facecolor="white")
    ax.set_facecolor("white")

    bp = ax.boxplot(
        grouped,
        positions=np.arange(1, len(order) + 1),
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": edge, "linewidth": 2.0},
        whiskerprops={"linewidth": 1.7},
        capprops={"linewidth": 1.7},
        boxprops={"linewidth": 1.8},
    )

    for i, name in enumerate(order):
        bp["boxes"][i].set_facecolor(fill[name])
        bp["boxes"][i].set_edgecolor(edge)
        bp["boxes"][i].set_alpha(0.92)
        for whisker in bp["whiskers"][2 * i : 2 * i + 2]:
            whisker.set_color(whisk[name])
        for cap in bp["caps"][2 * i : 2 * i + 2]:
            cap.set_color(whisk[name])

    rng = np.random.default_rng(42)
    for i, name in enumerate(order, start=1):
        vals = fit_df.loc[fit_df["model"] == name, "auc"].to_numpy()
        jitter = rng.normal(0.0, 0.085, size=vals.shape[0])
        ax.scatter(
            np.full(vals.shape[0], i, dtype=float) + jitter,
            vals,
            s=8,
            color="#5A5A5A",
            alpha=0.30,
            linewidths=0,
            zorder=1,
        )

    ax.set_title("Machine Learning model benchmark", fontsize=24, fontweight="bold", pad=18)
    ax.set_ylabel("Mean AUC", fontsize=20)
    ax.set_xlabel("")
    ax.set_xticks(np.arange(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=55, ha="right", rotation_mode="anchor", fontsize=12)
    ax.tick_params(axis="y", labelsize=16)
    ax.grid(True, axis="both", color="#D9D9D9", linewidth=1.1, alpha=0.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    y_min = float(fit_df["auc"].min())
    y_max = float(fit_df["auc"].max())
    ax.set_ylim(max(0.84, y_min - 0.02), min(1.005, y_max + 0.015))
    ax.set_xlim(0.45, len(order) + 0.55)

    fig.subplots_adjust(left=0.10, right=0.99, top=0.90, bottom=0.25)
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> int:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    classifier_meta_path = args.meta_path_classifier or os.path.join(
        args.run_dir, "tables", "metadata_strict_latest13k.csv"
    )
    if not os.path.exists(classifier_meta_path):
        raise RuntimeError(f"Missing classifier metadata CSV: {classifier_meta_path}")

    X, y, genes = load_classifier_matrix(
        adata_path=args.adata_path,
        classifier_meta_path=classifier_meta_path,
        gene_file=args.gene_file,
        run_dir=args.run_dir,
    )
    X_train, _, y_train, _ = train_test_split(
        X,
        y,
        test_size=0.20,
        random_state=args.seed,
        stratify=y,
    )

    fit_df = repeated_benchmark_df(
        X_train=X_train,
        y_train=y_train,
        seed=args.seed,
        cv_splits=args.cv_splits,
        cv_repeats=args.cv_repeats,
    )
    fit_df.to_csv(os.path.join(args.out_dir, "Fig6_panel_D_benchmark_points.csv"), index=False)
    summary = (
        fit_df.groupby("model")["auc"]
        .agg(["min", "median", "mean", "max", "std", "count"])
        .round(6)
        .reset_index()
    )
    summary.to_csv(os.path.join(args.out_dir, "Fig6_panel_D_benchmark_summary.csv"), index=False)

    out_png = os.path.join(args.out_dir, "Fig6_panel_D_only.png")
    out_pdf = os.path.join(args.out_dir, "Fig6_panel_D_only.pdf")
    plot_panel_d(fit_df, out_png, out_pdf, dpi=args.dpi)

    print(f"[OK] Saved: {out_png}")
    print(f"[OK] Saved: {out_pdf}")
    print(f"[OK] Saved: {os.path.join(args.out_dir, 'Fig6_panel_D_benchmark_points.csv')}")
    print(f"[OK] Saved: {os.path.join(args.out_dir, 'Fig6_panel_D_benchmark_summary.csv')}")
    print(f"[INFO] Genes used ({len(genes)}): {genes}")
    print(f"[INFO] Classifier cohort total: {len(y)}")
    print(f"[INFO] Training cohort used for D: {len(y_train)}")
    print(f"[INFO] CV setting: {args.cv_splits}-fold x {args.cv_repeats} repeats")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
