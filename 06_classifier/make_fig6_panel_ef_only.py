#!/usr/bin/env python3
"""
Render only Figure 6 panel E/F from an existing fig6 model cache.
"""

from __future__ import annotations

import argparse
import glob
import os
import pickle
from typing import Dict, List, Tuple

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["pdf.use14corefonts"] = False
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.sans-serif"] = ["Arial"]
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, precision_recall_curve, roc_curve


DEFAULT_RUN_DIR = (
    "results/reproduce_paper/runs/significant_noclobber_full"
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Render Figure 6 panel E/F only.")
    p.add_argument(
        "--run-dir",
        default=DEFAULT_RUN_DIR,
        help="Run directory containing fig6_model_cache_*.pkl.",
    )
    p.add_argument(
        "--cache-path",
        default="",
        help="Optional explicit fig6 model cache pickle.",
    )
    p.add_argument(
        "--out-dir",
        default="results/reproduce_paper/figures/fig6_panel_ef_only",
        help="Directory for outputs.",
    )
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def resolve_cache_path(run_dir: str, cache_path: str) -> str:
    if cache_path:
        if not os.path.exists(cache_path):
            raise RuntimeError(f"Missing cache path: {cache_path}")
        return cache_path
    hits = sorted(glob.glob(os.path.join(run_dir, "tables", "fig6_model_cache_*.pkl")))
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise RuntimeError(f"No fig6_model_cache_*.pkl found under {run_dir}/tables")
    raise RuntimeError(f"Multiple cache files found; pass --cache-path explicitly: {hits}")


def load_cache(path: str) -> dict:
    with open(path, "rb") as f:
        return pickle.load(f)


def ef_model_order() -> List[str]:
    return [
        "classif.log_reg",
        "classif.xgboost",
        "classif.lda",
        "classif.svm",
        "classif.naive_bayes",
        "classif.kknn",
        "classif.rpart",
        "classif.ranger",
    ]


def ef_model_palette() -> Dict[str, str]:
    return {
        "classif.log_reg": "#5B2A86",
        "classif.xgboost": "#4D3F99",
        "classif.lda": "#4965B2",
        "classif.svm": "#348DBE",
        "classif.naive_bayes": "#37AFC4",
        "classif.kknn": "#40C7C0",
        "classif.rpart": "#4ECE83",
        "classif.ranger": "#97DB46",
    }


def style_curve_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.grid(True, color="#D8D8D8", alpha=0.45, linewidth=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#BEBEBE")
    ax.spines["bottom"].set_color("#BEBEBE")
    ax.tick_params(labelsize=11, colors="#5C5C5C")


def plot_panel_e(cache: dict, out_png: str, out_pdf: str, dpi: int) -> pd.DataFrame:
    fig = plt.figure(figsize=(12.0, 4.9), facecolor="white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.42], wspace=0.28)
    ax_pr = fig.add_subplot(gs[0, 0])
    ax_roc = fig.add_subplot(gs[0, 1])
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_leg.axis("off")

    y_test = np.asarray(cache["y_test"], dtype=int)
    palette = ef_model_palette()
    rows = []

    for name in ef_model_order():
        rec = next(r for r in cache["fit_records"] if r["name"] == name)
        proba = np.asarray(rec["test_proba"], dtype=float)
        prec, recall, _ = precision_recall_curve(y_test, proba)
        fpr, tpr, _ = roc_curve(y_test, proba)
        roc_auc = float(auc(fpr, tpr))
        pr_auc = float(auc(recall, prec))
        rows.append({"model": name, "roc_auc": roc_auc, "pr_auc": pr_auc})
        ax_pr.plot(recall, prec, color=palette[name], lw=1.7, label=name)
        ax_roc.plot(fpr, tpr, color=palette[name], lw=1.7, label=name)

    style_curve_axes(ax_pr)
    style_curve_axes(ax_roc)
    ax_pr.axhline(0.5, color="#BFBFBF", lw=1.0, ls="--", zorder=0)
    ax_roc.plot([0, 1], [0, 1], color="#BFBFBF", lw=1.0, ls="--", zorder=0)

    ax_pr.set_xlim(0, 1.02)
    ax_pr.set_ylim(0, 1.04)
    ax_pr.set_xlabel("Recall", fontsize=16)
    ax_pr.set_ylabel("Precision", fontsize=16)
    ax_pr.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

    ax_roc.set_xlim(0, 1.02)
    ax_roc.set_ylim(0, 1.04)
    ax_roc.set_xlabel("1 - Specificity", fontsize=16)
    ax_roc.set_ylabel("Sensitivity", fontsize=16)
    ax_roc.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

    handles, labels = ax_roc.get_legend_handles_labels()
    ax_leg.legend(
        handles,
        labels,
        title="Learner",
        loc="center left",
        frameon=False,
        fontsize=11,
        title_fontsize=13,
        handlelength=1.4,
        handletextpad=0.6,
        borderaxespad=0.0,
    )

    fig.text(0.012, 0.985, "E", fontsize=22, fontweight="bold", ha="left", va="top")
    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.20)
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return pd.DataFrame(rows)


def style_small_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    for side in ax.spines.values():
        side.set_color("#AFAFAF")
        side.set_linewidth(0.9)
    ax.tick_params(labelsize=8, colors="#555555", length=2.5)


def plot_panel_f(cache: dict, out_png: str, out_pdf: str, dpi: int) -> pd.DataFrame:
    fig = plt.figure(figsize=(7.0, 6.6), facecolor="white")
    gs = fig.add_gridspec(2, 2, wspace=0.34, hspace=0.34)
    ax11 = fig.add_subplot(gs[0, 0])
    ax12 = fig.add_subplot(gs[0, 1])
    ax21 = fig.add_subplot(gs[1, 0])
    ax22 = fig.add_subplot(gs[1, 1])

    y_train = np.asarray(cache["y_train"], dtype=int)
    y_test = np.asarray(cache["y_test"], dtype=int)
    best_name = str(cache["best_name"])
    rec = next(r for r in cache["fit_records"] if r["name"] == best_name)
    train_proba = np.asarray(rec["train_proba"], dtype=float)
    test_proba = np.asarray(rec["test_proba"], dtype=float)

    fpr_tr, tpr_tr, _ = roc_curve(y_train, train_proba)
    fpr_te, tpr_te, _ = roc_curve(y_test, test_proba)
    prec_tr, rec_tr, _ = precision_recall_curve(y_train, train_proba)
    prec_te, rec_te, _ = precision_recall_curve(y_test, test_proba)
    auc_tr = float(auc(fpr_tr, tpr_tr))
    auc_te = float(auc(fpr_te, tpr_te))
    pr_tr = float(auc(rec_tr, prec_tr))
    pr_te = float(auc(rec_te, prec_te))
    line_color = "#5FC7CC"

    ax11.plot(fpr_tr, tpr_tr, color=line_color, lw=1.35)
    ax11.plot([0, 1], [0, 1], color="#BFBFBF", lw=0.9, ls="--")
    ax11.set_title(f"development set N={len(y_train)}, AUC = {auc_tr:.2f}", fontsize=8, pad=4)
    ax11.set_xlabel("False Positive Rate (FPR)", fontsize=9)
    ax11.set_ylabel("True Positive Rate (TPR)", fontsize=9)

    ax12.plot(rec_tr, prec_tr, color=line_color, lw=1.35)
    ax12.axhline(float(np.mean(y_train)), color="#BFBFBF", lw=0.9, ls="--")
    ax12.set_title(f"development set N={len(y_train)}, PRAUC = {pr_tr:.2f}", fontsize=8, pad=4)
    ax12.set_xlabel("Recall", fontsize=9)
    ax12.set_ylabel("Precision", fontsize=9)

    ax21.plot(fpr_te, tpr_te, color=line_color, lw=1.35)
    ax21.plot([0, 1], [0, 1], color="#BFBFBF", lw=0.9, ls="--")
    ax21.set_title(f"test set N={len(y_test)}, AUC = {auc_te:.2f}", fontsize=8, pad=4)
    ax21.set_xlabel("False Positive Rate (FPR)", fontsize=9)
    ax21.set_ylabel("True Positive Rate (TPR)", fontsize=9)

    ax22.plot(rec_te, prec_te, color=line_color, lw=1.35)
    ax22.axhline(float(np.mean(y_test)), color="#BFBFBF", lw=0.9, ls="--")
    ax22.set_title(f"test set N={len(y_test)}, PRAUC = {pr_te:.2f}", fontsize=8, pad=4)
    ax22.set_xlabel("Recall", fontsize=9)
    ax22.set_ylabel("Precision", fontsize=9)

    for ax in [ax11, ax12, ax21, ax22]:
        style_small_axes(ax)
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.02)

    fig.text(0.01, 0.985, "F", fontsize=22, fontweight="bold", ha="left", va="top")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.10)
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return pd.DataFrame(
        [
            {"split": "development", "metric": "ROC_AUC", "value": auc_tr},
            {"split": "development", "metric": "PR_AUC", "value": pr_tr},
            {"split": "test", "metric": "ROC_AUC", "value": auc_te},
            {"split": "test", "metric": "PR_AUC", "value": pr_te},
            {"split": "best_model", "metric": "name", "value": best_name},
        ]
    )


def main() -> int:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    cache_path = resolve_cache_path(args.run_dir, args.cache_path)
    cache = load_cache(cache_path)

    e_png = os.path.join(args.out_dir, "Fig6_panel_E_only.png")
    e_pdf = os.path.join(args.out_dir, "Fig6_panel_E_only.pdf")
    f_png = os.path.join(args.out_dir, "Fig6_panel_F_only.png")
    f_pdf = os.path.join(args.out_dir, "Fig6_panel_F_only.pdf")

    e_df = plot_panel_e(cache, e_png, e_pdf, args.dpi)
    f_df = plot_panel_f(cache, f_png, f_pdf, args.dpi)
    e_df.to_csv(os.path.join(args.out_dir, "Fig6_panel_E_metrics.csv"), index=False)
    f_df.to_csv(os.path.join(args.out_dir, "Fig6_panel_F_metrics.csv"), index=False)

    print(f"[OK] Saved: {e_png}")
    print(f"[OK] Saved: {e_pdf}")
    print(f"[OK] Saved: {f_png}")
    print(f"[OK] Saved: {f_pdf}")
    print(f"[OK] Saved: {os.path.join(args.out_dir, 'Fig6_panel_E_metrics.csv')}")
    print(f"[OK] Saved: {os.path.join(args.out_dir, 'Fig6_panel_F_metrics.csv')}")
    print(f"[INFO] Cache used: {cache_path}")
    print(f"[INFO] Best model: {cache['best_name']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
