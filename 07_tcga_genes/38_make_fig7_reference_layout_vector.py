#!/usr/bin/env python3
"""
Compose Figure 7 using manual-extracted reference coordinates (A1..D3).
Outputs fully editable vector files (SVG/PDF) plus optional PNG preview.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.metrics import roc_auc_score, roc_curve

# Keep text editable in vector exports.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build Figure 7 with reference panel coordinates (editable vector)."
    )
    p.add_argument("--run-dir", required=True, help="Run directory (reads tables/, writes figures/logs).")
    p.add_argument("--layout-meta", required=True, help="Path to metadata.json from panel extraction.")
    p.add_argument(
        "--expr-file",
        default="data/external/tcga_luad/expression.tsv.gz",
        help="Expression TSV(.gz), genes x samples.",
    )
    p.add_argument(
        "--surv-file",
        default="data/external/tcga_luad/survival.tsv.gz",
        help="Survival TSV(.gz), requires sample/OS/OS.time.",
    )
    p.add_argument(
        "--map-file",
        default="data/external/ensembl_to_symbol.csv",
        help="Ensembl->symbol mapping CSV with gene_id,gene_name.",
    )
    p.add_argument(
        "--genes-file",
        default="",
        help="Gene file with `gene` column. Default: <run-dir>/tables/final_consensus_genes.csv",
    )
    p.add_argument("--tumor-label", default="LUAD", help="Tumor group label in A/B panels.")
    p.add_argument("--cohort-label", default="TCGA-LUAD", help="Cohort label for KM titles.")
    p.add_argument("--xmax-years", type=float, default=10.0, help="KM x-axis max in years.")
    p.add_argument(
        "--out-prefix",
        default="Fig7_reference_layout_follow_fig5_vector",
        help="Output prefix under <run-dir>/figures and logs.",
    )
    p.add_argument("--png-dpi", type=int, default=300, help="Preview PNG DPI.")
    p.add_argument(
        "--canvas-width-in",
        type=float,
        default=13.5,
        help="Physical canvas width (inches). Height follows layout aspect ratio.",
    )
    p.add_argument(
        "--ai-optimize",
        action="store_true",
        help="Apply typography/spacing optimization for AI editing layout.",
    )
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


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "P = NA"
    if p < 1e-4:
        return f"P = {p:.2e}"
    return f"P = {p:.4f}"


def load_expression(expr_file: str, map_file: str) -> pd.DataFrame:
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


def prep_data(expr: pd.DataFrame, surv_file: str, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    genes_avail = [g for g in genes if g in expr.index]
    if not genes_avail:
        raise ValueError("None of requested genes found in expression matrix.")

    expr_sub = expr.loc[genes_avail].T.copy()
    expr_sub[genes_avail] = np.log2(expr_sub[genes_avail].clip(lower=0) + 1.0)
    expr_sub["sample"] = expr_sub.index.astype(str)
    expr_sub["sample_short"] = expr_sub["sample"].str[:15]
    expr_sub["patient"] = expr_sub["sample"].str[:12]
    expr_sub["type"] = expr_sub["sample"].apply(infer_sample_type)
    expr_use = expr_sub[expr_sub["type"].isin(["Normal", "Tumor"])].copy()

    surv = pd.read_csv(surv_file, sep="\t")
    surv["sample_short"] = surv["sample"].astype(str).str[:15]
    tumor_expr = expr_use[expr_use["type"] == "Tumor"][["sample_short"] + genes_avail].copy()
    merged = pd.merge(surv, tumor_expr, on="sample_short", how="inner")
    merged = merged.dropna(subset=["OS", "OS.time"]).copy()
    merged["OS"] = merged["OS"].astype(int)
    merged["OS.time"] = merged["OS.time"].astype(float)
    merged["time_year"] = merged["OS.time"] / 365.0
    return expr_use, merged, genes_avail


def style_axis(ax: plt.Axes, tick_size: float = 8.0) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=tick_size)


def panel_a(ax: plt.Axes, expr_use: pd.DataFrame, gene: str, tumor_label: str) -> float:
    normal = expr_use.loc[expr_use["type"] == "Normal", gene].astype(float)
    tumor = expr_use.loc[expr_use["type"] == "Tumor", gene].astype(float)
    if normal.empty or tumor.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        return np.nan

    data = [normal.values, tumor.values]
    parts = ax.violinplot(data, positions=[0, 1], widths=0.85, showextrema=False)
    for i, body in enumerate(parts["bodies"]):
        body.set_alpha(0.38)
        body.set_facecolor("#149C9A" if i == 0 else "#D87D4A")
        body.set_edgecolor("none")

    bp = ax.boxplot(data, positions=[0, 1], widths=0.15, patch_artist=True, showfliers=False)
    bp["boxes"][0].set(facecolor="#BFDCDC", edgecolor="black", alpha=0.9)
    bp["boxes"][1].set(facecolor="#F2B38E", edgecolor="black", alpha=0.9)
    for k in ["whiskers", "caps", "medians"]:
        for it in bp[k]:
            it.set(color="black", linewidth=0.9)

    rng = np.random.default_rng(42)
    ax.scatter(rng.normal(0, 0.05, len(normal)), normal, s=5, color="#108A88", linewidths=0)
    ax.scatter(rng.normal(1, 0.05, len(tumor)), tumor, s=5, color="#C81D25", linewidths=0)

    p = mannwhitneyu(normal, tumor, alternative="two-sided").pvalue
    ymin = min(float(normal.min()), float(tumor.min()))
    ymax = max(float(normal.max()), float(tumor.max()))
    span = ymax - ymin if ymax > ymin else 1.0
    ax.text(0.5, ymax + 0.075 * span, fmt_p(p), ha="center", va="bottom", fontsize=8.5)

    ax.set_title(gene, fontsize=12.5, fontweight="bold", pad=3, y=1.01)
    ax.set_xticks([0, 1], ["Normal", tumor_label])
    ax.set_ylabel(gene, fontsize=9.2, labelpad=1)
    ax.set_ylim(ymin - 0.08 * span, ymax + 0.18 * span)
    style_axis(ax, tick_size=8.5)
    return float(p)


def panel_b(ax: plt.Axes, expr_use: pd.DataFrame, gene: str, tumor_label: str) -> tuple[float, int]:
    pivot = (
        expr_use.pivot_table(index="patient", columns="type", values=gene, aggfunc="mean")
        .dropna()
        .sort_index()
    )
    if pivot.empty:
        ax.text(0.5, 0.5, "No paired data", ha="center", va="center")
        return np.nan, 0

    for _, row in pivot.iterrows():
        ax.plot([0, 1], [row["Normal"], row["Tumor"]], color="black", alpha=0.72, lw=0.9)
        ax.scatter([0, 1], [row["Normal"], row["Tumor"]], color="black", s=6, zorder=3)

    data = [pivot["Normal"].values, pivot["Tumor"].values]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.42, patch_artist=True, showfliers=False)
    bp["boxes"][0].set(facecolor="#149C9A", edgecolor="black", alpha=0.9)
    bp["boxes"][1].set(facecolor="#C81D25", edgecolor="black", alpha=0.9)
    for k in ["whiskers", "caps", "medians"]:
        for it in bp[k]:
            it.set(color="black", linewidth=0.9)

    try:
        p = wilcoxon(pivot["Normal"], pivot["Tumor"]).pvalue
    except ValueError:
        p = np.nan

    ymin = float(np.nanmin(pivot.values))
    ymax = float(np.nanmax(pivot.values))
    span = ymax - ymin if ymax > ymin else 1.0
    ax.text(0.5, ymax + 0.06 * span, fmt_p(p), ha="center", va="bottom", fontsize=8.5)

    ax.set_title(gene, fontsize=12.5, fontweight="bold", pad=3, y=1.01)
    ax.set_xticks([0, 1], ["Normal", tumor_label])
    ax.set_ylabel(f"{gene} expression", fontsize=9.2, labelpad=1)
    ax.set_ylim(ymin - 0.08 * span, ymax + 0.15 * span)
    style_axis(ax, tick_size=8.5)
    return float(p), int(len(pivot))


def panel_c1(ax: plt.Axes, expr_use: pd.DataFrame, genes: list[str]) -> dict[str, float]:
    y = (expr_use["type"] == "Tumor").astype(int).values
    palette = ["#273C75", "#C0392B", "#D4A62A", "#2E8B57", "#6AB7BE", "#8E44AD"]
    auc_map: dict[str, float] = {}
    for i, gene in enumerate(genes):
        score = expr_use[gene].astype(float).values
        if len(np.unique(y)) < 2 or len(np.unique(score)) < 2:
            auc_map[gene] = np.nan
            continue
        fpr, tpr, _ = roc_curve(y, score)
        auc = roc_auc_score(y, score)
        auc_map[gene] = float(auc)
        ax.plot(fpr, tpr, lw=2.2, color=palette[i % len(palette)], label=f"{gene}, AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], color="#BBBBBB", lw=1.0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("1 - Specificity", fontsize=11.5)
    ax.set_ylabel("Sensitivity", fontsize=11.5, labelpad=5)
    ax.tick_params(axis="both", labelsize=9)
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)
    style_axis(ax, tick_size=9.2)
    return auc_map


def at_risk(times: np.ndarray, cut_times: list[float]) -> list[int]:
    return [int(np.sum(times >= t)) for t in cut_times]


def km_block(
    fig: plt.Figure,
    rect: tuple[float, float, float, float],
    merged_tumor: pd.DataFrame,
    gene: str,
    cohort_label: str,
    xmax_years: float,
    show_left_labels: bool = True,
) -> tuple[float, int, int]:
    left, bottom, width, height = rect
    main_h = height * 0.76
    risk_h = height * 0.21
    gap = height * 0.03
    ax = fig.add_axes([left, bottom + risk_h + gap, width, main_h])
    tax = fig.add_axes([left, bottom, width, risk_h])

    df = merged_tumor[[gene, "OS", "time_year"]].dropna().copy()
    if df.empty:
        ax.text(0.5, 0.5, "No survival data", ha="center", va="center")
        tax.axis("off")
        return np.nan, 0, 0

    med = df[gene].median()
    df["Group"] = np.where(df[gene] > med, "High", "Low")
    high = df[df["Group"] == "High"]
    low = df[df["Group"] == "Low"]
    if high.empty or low.empty:
        ax.text(0.5, 0.5, "Empty group", ha="center", va="center")
        tax.axis("off")
        return np.nan, len(high), len(low)

    kmf = KaplanMeierFitter()
    kmf.fit(high["time_year"], high["OS"], label="High")
    kmf.plot_survival_function(
        ax=ax,
        ci_show=False,
        color="#C81D25",
        lw=1.8,
        show_censors=True,
        censor_styles={"marker": "|", "ms": 6},
    )
    kmf.fit(low["time_year"], low["OS"], label="Low")
    kmf.plot_survival_function(
        ax=ax,
        ci_show=False,
        color="#128A8A",
        lw=1.8,
        show_censors=True,
        censor_styles={"marker": "|", "ms": 6},
    )

    p = logrank_test(
        high["time_year"],
        low["time_year"],
        event_observed_A=high["OS"],
        event_observed_B=low["OS"],
    ).p_value

    cut_times = [0, 2.5, 5, 7.5, 10]

    ax.set_title(f"{cohort_label}: {gene}", fontsize=11.3, fontweight="bold", pad=2.5, y=1.01)
    ax.set_xlim(0, xmax_years)
    ax.set_ylim(0, 1.05)
    ax.set_xticks(cut_times)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelbottom=False)
    ax.set_ylabel("Survival probability", fontsize=8.4, fontweight="bold", labelpad=1)
    ax.legend(
        title=None,
        loc="upper right",
        frameon=False,
        fontsize=7.0,
        ncol=2,
        handlelength=1.1,
        columnspacing=0.6,
    )
    ax.text(0.04 * xmax_years, 0.11, f"p = {p:.3g}", fontsize=7.8)
    style_axis(ax, tick_size=7.8)

    high_risk = at_risk(high["time_year"].values, cut_times)
    low_risk = at_risk(low["time_year"].values, cut_times)

    tax.set_xlim(0, xmax_years)
    tax.set_ylim(0, 1)
    tax.set_yticks([])
    tax.set_xticks(cut_times)
    tax.tick_params(axis="x", labelsize=6.8, pad=1)
    tax.set_xlabel("Time(years)", fontsize=7.2, fontweight="bold", labelpad=2)
    for s in ["top", "right", "left"]:
        tax.spines[s].set_visible(False)

    if show_left_labels:
        tax.text(-0.055, 0.50, "mRNA", transform=tax.transAxes, rotation=90, va="center", ha="right", fontsize=7.4, fontweight="bold", clip_on=False)
        tax.text(-0.004, 0.72, "High", transform=tax.transAxes, ha="right", va="center", fontsize=6.9, clip_on=False)
        tax.text(-0.004, 0.30, "Low", transform=tax.transAxes, ha="right", va="center", fontsize=6.9, clip_on=False)
    for t, n in zip(cut_times, high_risk):
        tx = 0.35 if t == 0 else t
        tax.text(tx, 0.72, str(n), ha="center", va="center", fontsize=6.6)
    for t, n in zip(cut_times, low_risk):
        tx = 0.35 if t == 0 else t
        tax.text(tx, 0.30, str(n), ha="center", va="center", fontsize=6.6)

    return float(p), int(len(high)), int(len(low))


def to_rect(panel: dict, fig_w: float, fig_h: float) -> tuple[float, float, float, float]:
    x = float(panel["x"]) / fig_w
    y = float(panel["y"]) / fig_h
    w = float(panel["w"]) / fig_w
    h = float(panel["h"]) / fig_h
    bottom = 1.0 - y - h
    return x, bottom, w, h


def inset_rect(
    rect: tuple[float, float, float, float],
    left: float = 0.0,
    right: float = 0.0,
    bottom: float = 0.0,
    top: float = 0.0,
) -> tuple[float, float, float, float]:
    x, y, w, h = rect
    nx = x + w * left
    ny = y + h * bottom
    nw = w * max(0.05, 1.0 - left - right)
    nh = h * max(0.05, 1.0 - bottom - top)
    return nx, ny, nw, nh


def aligned_panel_rect(
    panel_map: dict[str, dict],
    panel_id: str,
    row_ids: list[str],
    fig_w: float,
    fig_h: float,
    align_row: bool = False,
) -> tuple[float, float, float, float]:
    p = dict(panel_map[panel_id])
    if align_row:
        row_top = min(float(panel_map[r]["y"]) for r in row_ids)
        row_bottom = max(float(panel_map[r]["y"]) + float(panel_map[r]["h"]) for r in row_ids)
        p["y"] = row_top
        p["h"] = row_bottom - row_top
    return to_rect(p, fig_w, fig_h)


def row_letter_pos(
    panel_map: dict[str, dict], panel_id: str, fig_w: float, fig_h: float
) -> tuple[float, float]:
    p = panel_map[panel_id]
    x = float(p["x"]) / fig_w
    y = float(p["y"]) / fig_h
    lx = max(0.003, x - 0.014)
    ly = min(0.995, 1.0 - y - 0.003)
    return lx, ly


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir).expanduser().resolve()
    fig_dir = run_dir / "figures"
    log_dir = run_dir / "logs"
    fig_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    layout_meta = Path(args.layout_meta).expanduser().resolve()
    if not layout_meta.exists():
        raise FileNotFoundError(f"layout metadata not found: {layout_meta}")
    with open(layout_meta, "r", encoding="utf-8") as f:
        meta = json.load(f)
    fig_w_px = float(meta["figure"]["w"])
    fig_h_px = float(meta["figure"]["h"])
    panel_map = {p["panel_id"]: p for p in meta["panels"]}

    required = [f"A{i}" for i in range(1, 6)] + [f"B{i}" for i in range(1, 6)] + ["C1", "C2", "C3", "D1", "D2", "D3"]
    missing = [k for k in required if k not in panel_map]
    if missing:
        raise ValueError(f"Missing panels in metadata: {missing}")

    if args.genes_file.strip():
        gene_file = Path(args.genes_file).expanduser().resolve()
    else:
        gene_file = run_dir / "tables" / "final_consensus_genes.csv"
    if not gene_file.exists():
        raise FileNotFoundError(f"gene file not found: {gene_file}")
    gdf = pd.read_csv(gene_file)
    if "gene" not in gdf.columns:
        raise ValueError(f"gene file missing `gene` column: {gene_file}")
    genes = [str(x).strip() for x in gdf["gene"].tolist() if str(x).strip()]
    if len(genes) != 5:
        raise ValueError(f"This reference layout expects 5 genes, got {len(genes)}: {genes}")

    expr = load_expression(args.expr_file, args.map_file)
    expr_use, merged_tumor, genes_avail = prep_data(expr, args.surv_file, genes)
    if len(genes_avail) != 5:
        raise ValueError(f"Expected 5 available genes, got {len(genes_avail)}: {genes_avail}")
    # Map genes to C/D blocks: C2,C3,D1,D2,D3
    km_gene_order = [genes_avail[0], genes_avail[1], genes_avail[2], genes_avail[3], genes_avail[4]]
    km_panel_order = ["C2", "C3", "D1", "D2", "D3"]
    row_ids = {
        "A": [f"A{i}" for i in range(1, 6)],
        "B": [f"B{i}" for i in range(1, 6)],
        "C": ["C1", "C2", "C3"],
        "D": ["D1", "D2", "D3"],
    }

    canvas_w = float(args.canvas_width_in)
    canvas_h = canvas_w * (fig_h_px / fig_w_px)
    fig = plt.figure(figsize=(canvas_w, canvas_h), dpi=100)
    fig.patch.set_facecolor("white")

    stats_a = []
    for i, g in enumerate(genes_avail, start=1):
        rect = aligned_panel_rect(
            panel_map, f"A{i}", row_ids["A"], fig_w_px, fig_h_px, align_row=args.ai_optimize
        )
        if args.ai_optimize:
            left = 0.14 if i == 1 else 0.09
            right = 0.06
            rect = inset_rect(rect, left=left, right=right, bottom=0.11, top=0.11)
        ax = fig.add_axes(rect)
        p = panel_a(ax, expr_use, g, args.tumor_label)
        stats_a.append({"gene": g, "panel": "A", "p_value": p})

    stats_b = []
    for i, g in enumerate(genes_avail, start=1):
        rect = aligned_panel_rect(
            panel_map, f"B{i}", row_ids["B"], fig_w_px, fig_h_px, align_row=args.ai_optimize
        )
        if args.ai_optimize:
            left = 0.14 if i == 1 else 0.09
            right = 0.06
            rect = inset_rect(rect, left=left, right=right, bottom=0.11, top=0.11)
        ax = fig.add_axes(rect)
        p, n = panel_b(ax, expr_use, g, args.tumor_label)
        stats_b.append({"gene": g, "panel": "B", "p_value": p, "n_pairs": n})

    rect_c1 = aligned_panel_rect(
        panel_map, "C1", row_ids["C"], fig_w_px, fig_h_px, align_row=args.ai_optimize
    )
    if args.ai_optimize:
        rect_c1 = inset_rect(rect_c1, left=0.105, right=0.02, bottom=0.08, top=0.06)
    ax_c1 = fig.add_axes(rect_c1)
    auc_map = panel_c1(ax_c1, expr_use, genes_avail)
    stats_c = [{"gene": g, "panel": "C", "auc": auc_map.get(g, np.nan)} for g in genes_avail]

    stats_d = []
    for panel_id, gene in zip(km_panel_order, km_gene_order):
        row_key = "C" if panel_id in ("C2", "C3") else "D"
        rect = aligned_panel_rect(
            panel_map, panel_id, row_ids[row_key], fig_w_px, fig_h_px, align_row=args.ai_optimize
        )
        if args.ai_optimize:
            left = 0.090 if panel_id in ("C2", "D1") else 0.070
            right = 0.03 if panel_id == "D3" else 0.022
            bottom = 0.145 if row_key == "C" else 0.17
            rect = inset_rect(rect, left=left, right=right, bottom=bottom, top=0.055)
        show_left_labels = panel_id in ("C2", "D1")
        p, nh, nl = km_block(
            fig,
            rect,
            merged_tumor,
            gene,
            args.cohort_label,
            args.xmax_years,
            show_left_labels=show_left_labels,
        )
        stats_d.append({"gene": gene, "panel": "D", "logrank_p": p, "n_high": nh, "n_low": nl})

    # Figure row tags from actual row anchors.
    letter_pos = {
        "A": row_letter_pos(panel_map, "A1", fig_w_px, fig_h_px),
        "B": row_letter_pos(panel_map, "B1", fig_w_px, fig_h_px),
        "C": row_letter_pos(panel_map, "C1", fig_w_px, fig_h_px),
        "D": row_letter_pos(panel_map, "D1", fig_w_px, fig_h_px),
    }
    if args.ai_optimize:
        # Anchor letters to aligned row tops (avoids row drift after insets).
        for k, ids in row_ids.items():
            top = min(float(panel_map[i]["y"]) for i in ids)
            x0 = float(panel_map[ids[0]]["x"]) / fig_w_px
            letter_pos[k] = (max(0.003, x0 - 0.012), min(0.995, 1.0 - top / fig_h_px - 0.001))
    for lab, (x, y) in letter_pos.items():
        fig.text(x, y, lab, fontsize=22 if args.ai_optimize else 28, fontweight="bold", va="top")

    out_png = fig_dir / f"{args.out_prefix}.png"
    out_pdf = fig_dir / f"{args.out_prefix}.pdf"
    out_svg = fig_dir / f"{args.out_prefix}.svg"
    fig.savefig(out_png, dpi=args.png_dpi, bbox_inches=None)
    fig.savefig(out_pdf, bbox_inches=None)
    fig.savefig(out_svg, bbox_inches=None)
    plt.close(fig)

    stats_df = pd.concat([pd.DataFrame(stats_a), pd.DataFrame(stats_b), pd.DataFrame(stats_c), pd.DataFrame(stats_d)], ignore_index=True, sort=False)
    stats_csv = log_dir / f"{args.out_prefix}_stats.csv"
    stats_df.to_csv(stats_csv, index=False)

    summary_txt = log_dir / f"{args.out_prefix}_summary.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write(f"run_dir={run_dir}\n")
        f.write(f"layout_meta={layout_meta}\n")
        f.write(f"figure_w_h_px={int(fig_w_px)}x{int(fig_h_px)}\n")
        f.write(f"genes_requested={','.join(genes)}\n")
        f.write(f"genes_used={','.join(genes_avail)}\n")
        f.write(f"expr_samples_total={expr.shape[1]}\n")
        f.write(f"expr_samples_normal_tumor={len(expr_use)}\n")
        f.write(f"tumor_survival_rows={len(merged_tumor)}\n")
        f.write(f"out_png={out_png}\n")
        f.write(f"out_pdf={out_pdf}\n")
        f.write(f"out_svg={out_svg}\n")
        f.write(f"stats_csv={stats_csv}\n")

    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_svg}")
    print(f"Saved stats: {stats_csv}")
    print(f"Saved summary: {summary_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
