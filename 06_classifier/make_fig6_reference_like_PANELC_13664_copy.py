#!/usr/bin/env python3
"""
Build a Figure 6 style composite (A-K) from local project data.

Outputs:
  - results/reproduce_paper/figures/Fig6_combined_reference_like.png
  - results/reproduce_paper/figures/Fig6_combined_reference_like.pdf
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import anndata as ad
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
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import gaussian_kde, spearmanr
from scipy.spatial import ConvexHull, Delaunay, cKDTree

try:
    import seaborn as sns
except ImportError:
    sns = None

try:
    import shap
except ImportError:
    shap = None

try:
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
    from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        auc,
        confusion_matrix,
        precision_recall_curve,
        roc_auc_score,
        roc_curve,
    )
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from sklearn.naive_bayes import GaussianNB
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import MinMaxScaler, StandardScaler
    from sklearn.svm import SVC
    from sklearn.tree import DecisionTreeClassifier
except ImportError:
    LinearDiscriminantAnalysis = None
    HistGradientBoostingClassifier = None
    RandomForestClassifier = None
    LogisticRegression = None
    auc = None
    confusion_matrix = None
    precision_recall_curve = None
    roc_auc_score = None
    roc_curve = None
    StratifiedKFold = None
    train_test_split = None
    GaussianNB = None
    KNeighborsClassifier = None
    make_pipeline = None
    MinMaxScaler = None
    StandardScaler = None
    SVC = None
    DecisionTreeClassifier = None


DEFAULT_RUN_DIR = (
    "results/reproduce_paper/runs/paper_main_epiref_latest13k_originalflow_20260214_125353"
)
DEFAULT_UMAP_CSV = "results/reproduce_paper/r_data/seurat_umap_coords.csv"
DEFAULT_MALIGNANT_UMAP_CSV = "results/reproduce_paper/r_data/malignant_strict_harmony_umap_coords.csv"
DEFAULT_SEURAT_RDS = "data/processed/step1_seurat_harmony.rds"
FALLBACK_GENES = ["APMAP", "CDKN2A", "TSPAN6", "MAPRE1", "GLB1"]


@dataclass
class FitResult:
    name: str
    model: object
    cv_auc: List[float]
    test_proba: np.ndarray
    train_proba: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Figure 6 reference-like composite.")
    parser.add_argument(
        "--run-dir",
        default=DEFAULT_RUN_DIR,
        help="Run directory containing tables/final_consensus_genes.csv.",
    )
    parser.add_argument(
        "--gene-file",
        default="",
        help="Optional explicit gene file (CSV with column 'gene'). Overrides --run-dir table.",
    )
    parser.add_argument(
        "--adata-path",
        default="data/processed/discovery_malignant.h5ad",
        help="Path to atlas-scale h5ad used for panel A expression/coordinates.",
    )
    parser.add_argument(
        "--meta-path",
        default="results/reproduce_paper/r_data/metadata_strict.csv",
        help="Path to strict metadata CSV used to define malignant subset and labels.",
    )
    parser.add_argument(
        "--meta-path-classifier",
        default="",
        help=(
            "Optional strict metadata CSV used only for D-K classifier panels. "
            "Defaults to <run-dir>/tables/metadata_strict_latest13k.csv."
        ),
    )
    parser.add_argument(
        "--meta-path-panel-b",
        default="",
        help=(
            "Optional metadata CSV used only for panel B malignant-cell selection. "
            "Defaults to runfinal inferCNV export if it can be resolved from --run-dir; "
            "otherwise falls back to --meta-path."
        ),
    )
    parser.add_argument(
        "--umap-csv",
        default=DEFAULT_UMAP_CSV,
        help="Optional CSV with columns cell_id,UMAP_1,UMAP_2 (from Fig1 Seurat object).",
    )
    parser.add_argument(
        "--umap-csv-malignant",
        default=DEFAULT_MALIGNANT_UMAP_CSV,
        help=(
            "Optional malignant-subset UMAP CSV (cell_id,UMAP_1,UMAP_2) for panel B. "
            "If missing, panel B is recomputed from Seurat Harmony embeddings."
        ),
    )
    parser.add_argument(
        "--seurat-rds",
        default=DEFAULT_SEURAT_RDS,
        help="Seurat RDS with Harmony embeddings used to rebuild malignant-only UMAP for panel B.",
    )
    parser.add_argument(
        "--out-dir",
        default="results/reproduce_paper/figures",
        help="Directory for output figures.",
    )
    parser.add_argument(
        "--layout-meta",
        default="",
        help=(
            "Optional metadata.json from manual panel extraction. "
            "If set, the full Figure 6 layout uses extracted panel coordinates "
            "(e.g. A1-A5, B1-B5, C-K) instead of the internal GridSpec guess."
        ),
    )
    parser.add_argument(
        "--fig-width",
        type=float,
        default=24.0,
        help="Composite figure width in inches.",
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=17.5,
        help="Composite figure height in inches.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Raster DPI for PNG export.",
    )
    parser.add_argument(
        "--analysis-cache",
        default="",
        help=(
            "Optional pickle cache for D-K model fits and SHAP results. "
            "If omitted, a run-specific cache path is created under <run-dir>/tables/."
        ),
    )
    parser.add_argument(
        "--refresh-analysis-cache",
        action="store_true",
        help="Ignore any existing analysis cache and recompute model/SHAP results.",
    )
    parser.add_argument(
        "--preview-only-abc",
        action="store_true",
        help="Render only panels A/B/C for fast layout debugging without model/SHAP panels.",
    )
    parser.add_argument(
        "--preview-only-def",
        action="store_true",
        help="Render only panels D/E/F using cached or freshly computed model results.",
    )
    parser.add_argument(
        "--compact",
        action="store_true",
        help="Use tighter layout spacing to make the full A-K panel more compact.",
    )
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def to_dense(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def read_metadata_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "cell_id" in df.columns:
        df = df.set_index("cell_id")
    elif df.columns.size > 0:
        first_col = str(df.columns[0])
        if first_col.startswith("Unnamed:") or first_col in {"index", ""}:
            df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    return df


def infer_panel_b_meta_path(run_dir: str) -> str:
    direct = os.path.join(run_dir, "tables", "cnv_paper_exact_cell_metadata_runfinal_ref_mean_3sd.csv")
    if os.path.exists(direct):
        return direct

    prov = os.path.join(run_dir, "tables", "run_input_provenance.txt")
    if not os.path.exists(prov):
        return ""

    input_rds = ""
    with open(prov, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("Input Seurat RDS:"):
                input_rds = line.split(":", 1)[1].strip()
                break
    if not input_rds:
        return ""

    upstream_run = os.path.dirname(os.path.dirname(input_rds))
    cand = os.path.join(upstream_run, "tables", "cnv_paper_exact_cell_metadata_runfinal_ref_mean_3sd.csv")
    return cand if os.path.exists(cand) else ""


def panel_b_umap_cache_path(panel_b_meta_path: str, explicit_umap_csv: str) -> str:
    if explicit_umap_csv != DEFAULT_MALIGNANT_UMAP_CSV:
        return explicit_umap_csv
    if not panel_b_meta_path:
        return explicit_umap_csv
    stem = os.path.splitext(os.path.basename(panel_b_meta_path))[0]
    return os.path.join("results", "reproduce_paper", "r_data", f"{stem}_harmony_umap_coords.csv")


def load_genes(gene_file: str, adata: ad.AnnData) -> List[str]:
    if not gene_file or not os.path.exists(gene_file):
        return FALLBACK_GENES.copy()
    df = pd.read_csv(gene_file)
    if "gene" not in df.columns:
        return FALLBACK_GENES.copy()
    genes = [str(g).strip() for g in df["gene"].tolist() if str(g).strip()]
    if not genes:
        return FALLBACK_GENES.copy()

    # Keep genes that exist in expression matrix. Preserve order.
    adata_genes = set(map(str, adata.raw.var_names if adata.raw is not None else adata.var_names))
    filtered = [g for g in genes if g in adata_genes]
    return filtered if filtered else FALLBACK_GENES.copy()


def default_analysis_cache_path(run_dir: str, classifier_meta_path: str, genes: List[str], seed: int) -> str:
    key = hashlib.md5(
        "|".join(
            [
                os.path.abspath(run_dir),
                os.path.abspath(classifier_meta_path),
                ",".join(genes),
                str(seed),
            ]
        ).encode("utf-8")
    ).hexdigest()[:12]
    return os.path.join(run_dir, "tables", f"fig6_model_cache_{key}.pkl")


def save_analysis_cache(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_analysis_cache(path: str) -> dict:
    with open(path, "rb") as f:
        return pickle.load(f)


def weighted_density(
    ax,
    xy: np.ndarray,
    w: np.ndarray,
    title: str,
    xlim: Tuple[float, float] | None = None,
    ylim: Tuple[float, float] | None = None,
    outline_xy: List[np.ndarray] | np.ndarray | None = None,
) -> None:
    if xy.shape[0] < 30:
        ax.text(0.5, 0.5, "Too few points", ha="center", va="center")
        ax.set_title(title, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        return

    x = xy[:, 0]
    y = xy[:, 1]
    if xlim is None:
        x_min, x_max = np.percentile(x, [1, 99])
        x_pad = (x_max - x_min) * 0.08
        x0, x1 = x_min - x_pad, x_max + x_pad
    else:
        x0, x1 = float(xlim[0]), float(xlim[1])

    if ylim is None:
        y_min, y_max = np.percentile(y, [1, 99])
        y_pad = (y_max - y_min) * 0.08
        y0, y1 = y_min - y_pad, y_max + y_pad
    else:
        y0, y1 = float(ylim[0]), float(ylim[1])

    xx, yy = np.mgrid[
        x0:x1:160j,
        y0:y1:160j,
    ]
    grid = np.vstack([xx.ravel(), yy.ravel()])

    w_pos = np.maximum(w, 0.0)
    if float(np.sum(w_pos)) <= 0:
        w_pos = np.ones_like(w_pos)

    kde_w = gaussian_kde(np.vstack([x, y]), weights=w_pos)
    zz = kde_w(grid).reshape(xx.shape)
    if outline_xy is not None:
        if isinstance(outline_xy, np.ndarray):
            outline_iter = [outline_xy]
        else:
            outline_iter = outline_xy
        keep_mask = np.zeros(xx.shape, dtype=bool)
        valid_polys = []
        for poly in outline_iter:
            if poly is None or poly.shape[0] < 3:
                continue
            valid_polys.append(poly)
            path = Path(poly[:-1] if np.allclose(poly[0], poly[-1]) else poly)
            keep_mask |= path.contains_points(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)
        if valid_polys:
            zz = np.ma.array(zz, mask=~keep_mask)
            outline_iter = valid_polys
        else:
            outline_xy = None
    m = ax.contourf(xx, yy, zz, levels=18, cmap="magma")
    cbar = plt.colorbar(m, ax=ax, fraction=0.040, pad=0.018)
    cbar.ax.tick_params(labelsize=7)

    if outline_xy is not None:
        plotted = False
        for poly in outline_iter:
            if poly is None or poly.shape[0] < 3:
                continue
            ax.plot(
                poly[:, 0],
                poly[:, 1],
                color="white",
                lw=1.6,
                ls="--",
            )
            plotted = True
        if not plotted:
            outline_xy = None
    if outline_xy is None:
        kde_shape = gaussian_kde(np.vstack([x, y]))
        zz_shape = kde_shape(grid).reshape(xx.shape)
        level = float(np.quantile(zz_shape[np.isfinite(zz_shape)], 0.80))
        ax.contour(
            xx,
            yy,
            zz_shape,
            levels=[level],
            colors="white",
            linewidths=1.6,
            linestyles="--",
        )
    ax.set_facecolor("black")
    ax.set_title(title, fontsize=9, color="white", pad=2, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)


def model_zoo(seed: int) -> Dict[str, object]:
    return {
        "classif.log_reg": make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=5000, random_state=seed),
        ),
        "classif.xgboost": HistGradientBoostingClassifier(
            max_depth=4,
            learning_rate=0.05,
            max_iter=300,
            random_state=seed,
        ),
        "classif.lda": make_pipeline(StandardScaler(), LinearDiscriminantAnalysis()),
        "classif.svm": make_pipeline(StandardScaler(), SVC(probability=True, random_state=seed)),
        "classif.naive_bayes": GaussianNB(),
        "classif.kknn": make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=31)),
        "classif.rpart": DecisionTreeClassifier(max_depth=4, random_state=seed),
        "classif.ranger": RandomForestClassifier(
            n_estimators=400, random_state=seed, n_jobs=-1
        ),
    }


def net_benefit(y_true: np.ndarray, y_prob: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob).astype(float)
    n = float(y_true.size)
    out = []
    for t in thresholds:
        pred = (y_prob >= t).astype(int)
        tp = float(np.sum((pred == 1) & (y_true == 1)))
        fp = float(np.sum((pred == 1) & (y_true == 0)))
        out.append(tp / n - fp / n * (t / (1.0 - t)))
    return np.asarray(out)


def format_conf_mat_reference(ax, y_true: np.ndarray, y_pred: np.ndarray, title: str) -> None:
    # Match the paper-style 3x3 confusion matrix with row/column totals.
    cm2 = confusion_matrix(y_true, y_pred, labels=[1, 0]).astype(int)
    total = int(np.sum(cm2))
    # confusion_matrix returns rows=true, cols=pred.
    # The paper panel uses rows=prediction, cols=truth.
    tp, fn = int(cm2[0, 0]), int(cm2[0, 1])
    fp, tn = int(cm2[1, 0]), int(cm2[1, 1])
    row1, row0 = tp + fp, fn + tn
    col1, col0 = tp + fn, fp + tn
    mat = np.array(
        [
            [tp, fp, row1],
            [fn, tn, row0],
            [col1, col0, total],
        ],
        dtype=float,
    )
    pct = (mat / float(total) * 100.0) if total > 0 else np.zeros_like(mat)

    def scale01(v: float, lo: float, hi: float) -> float:
        if hi - lo <= 1e-9:
            return 1.0
        return float(np.clip((v - lo) / (hi - lo), 0.0, 1.0))

    def blend_hex(c0: str, c1: str, t: float) -> np.ndarray:
        a = np.asarray(mpl.colors.to_rgb(c0), dtype=float)
        b = np.asarray(mpl.colors.to_rgb(c1), dtype=float)
        return a * (1.0 - t) + b * t

    # Per-role color ramps to reproduce the paper look.
    diag_lo, diag_hi = "#BDE4BA", "#3FAE57"
    off_lo, off_hi = "#EEF5EA", "#D3E4D1"
    total_lo, total_hi = "#F3E2D1", "#F46810"
    grand = "#E5E5E5"

    diag_vals = np.asarray([pct[0, 0], pct[1, 1]], dtype=float)
    off_vals = np.asarray([pct[0, 1], pct[1, 0]], dtype=float)
    total_vals = np.asarray([pct[0, 2], pct[1, 2], pct[2, 0], pct[2, 1]], dtype=float)

    diag_min, diag_max = float(np.min(diag_vals)), float(np.max(diag_vals))
    off_min, off_max = float(np.min(off_vals)), float(np.max(off_vals))
    total_min, total_max = float(np.min(total_vals)), float(np.max(total_vals))

    rgb = np.zeros((3, 3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            if i == 2 and j == 2:
                rgb[i, j, :] = np.asarray(mpl.colors.to_rgb(grand), dtype=float)
            elif i == 2 or j == 2:
                t = scale01(float(pct[i, j]), total_min, total_max)
                rgb[i, j, :] = blend_hex(total_lo, total_hi, t)
            elif i == j:
                t = scale01(float(pct[i, j]), diag_min, diag_max)
                rgb[i, j, :] = blend_hex(diag_lo, diag_hi, t)
            else:
                t = scale01(float(pct[i, j]), off_min, off_max)
                rgb[i, j, :] = blend_hex(off_lo, off_hi, t)

    ax.imshow(rgb, interpolation="nearest", origin="upper", aspect="equal")
    ax.set_xticks(np.arange(-0.5, 3.0, 1.0), minor=True)
    ax.set_yticks(np.arange(-0.5, 3.0, 1.0), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.grid(False)
    ax.tick_params(which="minor", bottom=False, left=False)

    def fmt_pct(v: float) -> str:
        if abs(v - round(v)) < 0.05:
            return f"{int(round(v))}%"
        return f"{v:.1f}%"

    for i in range(3):
        for j in range(3):
            if i == 2 and j == 2:
                ax.text(
                    float(j),
                    float(i),
                    f"{int(mat[i, j])}",
                    ha="center",
                    va="center",
                    fontsize=11,
                    color="black",
                )
                continue
            ax.text(
                float(j),
                float(i) - 0.06,
                fmt_pct(float(pct[i, j])),
                ha="center",
                va="center",
                fontsize=12 if (i < 2 and j < 2) else 11,
                color="black",
            )
            ax.text(
                float(j),
                float(i) + 0.12,
                f"{int(mat[i, j])}",
                ha="center",
                va="center",
                fontsize=8 if (i < 2 and j < 2) else 7,
                color="black",
            )

    # Tiny row/column rates around 2x2 core, as in the reference panel.
    core = mat[:2, :2]
    row_sums = np.sum(core, axis=1)
    col_sums = np.sum(core, axis=0)
    for i in range(2):
        for j in range(2):
            row_rate = 100.0 * core[i, j] / row_sums[i] if row_sums[i] > 0 else 0.0
            col_rate = 100.0 * core[i, j] / col_sums[j] if col_sums[j] > 0 else 0.0
            ax.text(
                float(j) + 0.42,
                float(i),
                f"{row_rate:.1f}%",
                rotation=90,
                ha="center",
                va="center",
                fontsize=6,
                color="#333333",
                alpha=0.85,
            )
            ax.text(
                float(j),
                float(i) + 0.38,
                f"{col_rate:.1f}%",
                ha="center",
                va="center",
                fontsize=6,
                color="#333333",
                alpha=0.85,
            )

    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(2.5, -0.5)
    ax.set_title(title, fontsize=12, pad=10)
    ax.set_xlabel("Truth (observation)", fontsize=9)
    ax.set_ylabel("Prediction (model output)", fontsize=9)
    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 1, 2])
    ax.set_xticklabels(["1", "0", "Total"], fontsize=8)
    ax.set_yticklabels(["1", "0", "Total"], rotation=90, fontsize=8)
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.tick_params(axis="both", length=0, colors="#555555")
    for spine in ax.spines.values():
        spine.set_visible(False)


def figure6_model_order(fit_names: List[str]) -> List[str]:
    preferred = [
        "classif.kknn",
        "classif.lda",
        "classif.log_reg",
        "classif.naive_bayes",
        "classif.ranger",
        "classif.rpart",
        "classif.svm",
        "classif.xgboost",
    ]
    rest = [n for n in fit_names if n not in preferred]
    return [n for n in preferred if n in fit_names] + rest


def figure6_model_palette(model_order: List[str]) -> Dict[str, tuple]:
    ref_colors = {
        "classif.kknn": "#6C8CC4",
        "classif.lda": "#E2B84D",
        "classif.log_reg": "#8B78B8",
        "classif.naive_bayes": "#D98AA8",
        "classif.ranger": "#7ED37F",
        "classif.rpart": "#F0A24A",
        "classif.svm": "#8FD5E8",
        "classif.xgboost": "#D782D0",
    }
    return {name: ref_colors.get(name, "#6A9FB5") for name in model_order}


def add_panel_letter(ax, letter: str) -> None:
    ax.text(
        -0.04,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def figure6_panel_g_form(train: bool) -> Dict[str, str]:
    return {
        "title": "Confusion matrix for the train set" if train else "Confusion matrix for the test set",
        "xlabel": "Truth (observation)",
        "ylabel": "Prediction (model output)",
    }


def figure6_panel_h_form(n_test: int) -> Dict[str, str]:
    return {
        "title": "Decision curve analysis",
        "annotation": f"number of cells in the test set N={n_test}",
        "xlabel": "threshold",
        "ylabel": "value",
    }


def figure6_panel_h_line_styles() -> List[Dict[str, object]]:
    return [
        {
            "label": "proposed model",
            "color": "#FF5A5A",
            "lw": 1.05,
            "dashes": (2.6, 1.4),
            "marker": "o",
            "markersize": 2.0,
            "markerfacecolor": "#FF5A5A",
            "markeredgewidth": 0.0,
            "markevery": 3,
        },
        {
            "label": "classify all as 1",
            "color": "#CFCFCF",
            "lw": 0.95,
            "dashes": (2.2, 1.5),
            "marker": "o",
            "markersize": 1.6,
            "markerfacecolor": "#D9D9D9",
            "markeredgewidth": 0.0,
            "markevery": 3,
        },
        {
            "label": "classify all as 0",
            "color": "#4F4F4F",
            "lw": 0.95,
            "dashes": (2.2, 1.5),
            "marker": "o",
            "markersize": 1.6,
            "markerfacecolor": "#4F4F4F",
            "markeredgewidth": 0.0,
            "markevery": 4,
        },
        {
            "label": "random guess",
            "color": "#53F05B",
            "lw": 0.95,
            "dashes": (2.4, 1.4),
            "marker": "o",
            "markersize": 1.8,
            "markerfacecolor": "#53F05B",
            "markeredgewidth": 0.0,
            "markevery": 3,
        },
    ]


def figure6_panel_h_legend_handles() -> List[mpl.lines.Line2D]:
    return [
        mpl.lines.Line2D(
            [],
            [],
            color=style["color"],
            lw=style["lw"],
            linestyle=(0, style["dashes"]),
            marker=style["marker"],
            markersize=style["markersize"],
            markerfacecolor=style["markerfacecolor"],
            markeredgewidth=style["markeredgewidth"],
            label=style["label"],
        )
        for style in figure6_panel_h_line_styles()
    ]


def figure6_panel_i_form() -> Dict[str, str]:
    return {
        "title": "SHAP Feature Importance Analysis",
        "xlabel": "SHAP Value",
        "ylabel": "",
        "legend_title": "correct_prediction",
        "colorbar_label": "Scaled Feature Value",
    }


def plot_panel_h_reference(ax, y_test: np.ndarray, test_proba: np.ndarray) -> None:
    form = figure6_panel_h_form(len(y_test))
    styles = figure6_panel_h_line_styles()
    th = np.linspace(0.01, 0.99, 181)
    nb_model = net_benefit(y_test, test_proba, th)
    prevalence = float(np.mean(y_test))
    nb_all = prevalence - (1.0 - prevalence) * (th / (1.0 - th))

    # Match the paper's auxiliary baselines:
    # - classify all as 0: horizontal zero line
    # - random guess: declining green reference curve
    nb_random = prevalence * 0.5 - (1.0 - prevalence) * 0.5 * (th / (1.0 - th))
    nb_all_plot = np.where(nb_all >= 0.0, nb_all, np.nan)
    nb_random_plot = np.where(nb_random >= 0.0, nb_random, np.nan)
    nb_none = np.zeros_like(th)

    ax.plot(
        th,
        nb_model,
        color=styles[0]["color"],
        lw=styles[0]["lw"],
        linestyle=(0, styles[0]["dashes"]),
        dash_capstyle="round",
        marker=styles[0]["marker"],
        markersize=styles[0]["markersize"],
        markerfacecolor=styles[0]["markerfacecolor"],
        markeredgewidth=styles[0]["markeredgewidth"],
        markevery=styles[0]["markevery"],
        label=styles[0]["label"],
    )
    ax.plot(
        th,
        nb_all_plot,
        color=styles[1]["color"],
        lw=styles[1]["lw"],
        linestyle=(0, styles[1]["dashes"]),
        dash_capstyle="round",
        marker=styles[1]["marker"],
        markersize=styles[1]["markersize"],
        markerfacecolor=styles[1]["markerfacecolor"],
        markeredgewidth=styles[1]["markeredgewidth"],
        markevery=styles[1]["markevery"],
        label=styles[1]["label"],
    )
    ax.plot(
        th,
        nb_none,
        color=styles[2]["color"],
        lw=styles[2]["lw"],
        linestyle=(0, styles[2]["dashes"]),
        dash_capstyle="round",
        marker=styles[2]["marker"],
        markersize=styles[2]["markersize"],
        markerfacecolor=styles[2]["markerfacecolor"],
        markeredgewidth=styles[2]["markeredgewidth"],
        markevery=styles[2]["markevery"],
        label=styles[2]["label"],
    )
    ax.plot(
        th,
        nb_random_plot,
        color=styles[3]["color"],
        lw=styles[3]["lw"],
        linestyle=(0, styles[3]["dashes"]),
        dash_capstyle="round",
        marker=styles[3]["marker"],
        markersize=styles[3]["markersize"],
        markerfacecolor=styles[3]["markerfacecolor"],
        markeredgewidth=styles[3]["markeredgewidth"],
        markevery=styles[3]["markevery"],
        label=styles[3]["label"],
    )
    ax.axvline(0.5, color="#2C5BFF", linestyle=(0, (3.2, 1.6)), lw=1.25)

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.01, max(0.52, float(np.nanmax(nb_model)) + 0.01))
    ax.set_title("")
    ax.text(
        0.0,
        1.06,
        form["title"],
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11,
        color="black",
        clip_on=False,
    )
    ax.text(
        0.0,
        1.00,
        form["annotation"],
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="black",
        clip_on=False,
    )
    ax.set_xlabel(form["xlabel"], fontsize=9)
    ax.set_ylabel(form["ylabel"], fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(
        frameon=False,
        fontsize=7,
        loc="center left",
        bbox_to_anchor=(1.03, 0.40),
        handlelength=2.2,
        borderaxespad=0.0,
    )
    ax.grid(False)
    ax.set_box_aspect(0.95)


def plot_panel_i_reference(
    ax,
    shap_vals: np.ndarray,
    feature_vals: np.ndarray,
    genes: List[str],
    correct_mask: np.ndarray,
) -> None:
    form = figure6_panel_i_form()
    shap_vals = np.asarray(shap_vals, dtype=float)
    feature_vals = np.asarray(feature_vals, dtype=float)
    correct_mask = np.asarray(correct_mask).astype(bool)
    n_genes = shap_vals.shape[1]

    mean_abs = np.mean(np.abs(shap_vals), axis=0)
    order = np.argsort(mean_abs)  # ascending, most important at bottom like reference

    x_max = float(np.quantile(np.abs(shap_vals), 0.992))
    x_max = float(np.clip(x_max * 1.08, 0.35, 0.45))
    x_lim = (-x_max, x_max)

    # Light gray violin background to match paper look.
    datasets = [shap_vals[:, idx] for idx in order]
    v = ax.violinplot(
        datasets,
        positions=np.arange(n_genes),
        vert=False,
        widths=0.52,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in v["bodies"]:
        body.set_facecolor("#D3D3D3")
        body.set_edgecolor("#BDBDBD")
        body.set_alpha(0.22)
        body.set_linewidth(0.45)

    cmap = plt.get_cmap("rainbow")
    norm = mpl.colors.PowerNorm(gamma=1.35, vmin=0.0, vmax=1.0)
    rng = np.random.default_rng(42)
    for pos, idx in enumerate(order):
        x = shap_vals[:, idx]
        fv = feature_vals[:, idx]
        lo, hi = float(np.min(fv)), float(np.max(fv))
        if hi - lo <= 1e-12:
            c = np.full_like(fv, 0.5, dtype=float)
        else:
            c = np.clip((fv - lo) / (hi - lo), 0.0, 1.0)
        y = pos + rng.normal(0.0, 0.085, size=x.shape[0])

        ax.scatter(
            x[correct_mask],
            y[correct_mask],
            c=c[correct_mask],
            cmap=cmap,
            norm=norm,
            s=7,
            marker="o",
            alpha=0.68,
            linewidths=0.0,
        )
        ax.scatter(
            x[~correct_mask],
            y[~correct_mask],
            c=c[~correct_mask],
            cmap=cmap,
            norm=norm,
            s=7,
            marker="x",
            alpha=0.56,
            linewidths=0.42,
        )

        ax.text(
            x_lim[0] + 0.01 * (x_lim[1] - x_lim[0]),
            pos,
            f"{mean_abs[idx]:.3f}",
            ha="left",
            va="center",
            fontsize=7,
            color="#333333",
        )

    ax.axvline(0.0, color="#8A8A8A", lw=0.8)
    ax.set_xlim(x_lim)
    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels([genes[idx] for idx in order], fontsize=8)
    ax.set_xlabel(form["xlabel"], fontsize=9)
    ax.set_ylabel(form["ylabel"], fontsize=9)
    ax.set_title(form["title"], fontsize=11, pad=6)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="x", alpha=0.20, color="#BEBEBE")
    ax.grid(axis="y", visible=False)

    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(
        ax,
        width="4.5%",
        height="70%",
        loc="lower left",
        bbox_to_anchor=(1.28, 0.12, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0.0,
    )
    cb = plt.colorbar(sm, cax=cax)
    cb.set_label(form["colorbar_label"], fontsize=8, labelpad=4)
    cb.set_ticks([0, 1])
    cb.set_ticklabels(["Low", "High"])
    cb.ax.tick_params(labelsize=8)

    handles = [
        mpl.lines.Line2D([], [], color="black", marker="x", linestyle="None", markersize=4, label="FALSE"),
        mpl.lines.Line2D([], [], color="black", marker="o", linestyle="None", markersize=4, label="TRUE"),
    ]
    ax.legend(
        handles=handles,
        title=form["legend_title"],
        frameon=False,
        fontsize=7,
        title_fontsize=8,
        loc="center left",
        bbox_to_anchor=(1.02, 0.63),
    )
    ax.set_box_aspect(0.95)


def linear_fit_ci_band(x: np.ndarray, y: np.ndarray, x_grid: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_grid = np.asarray(x_grid, dtype=float)
    n = x.shape[0]
    if n < 3 or float(np.std(x)) < 1e-10:
        mid = np.full_like(x_grid, float(np.mean(y)))
        return mid, mid, mid

    x_bar = float(np.mean(x))
    y_bar = float(np.mean(y))
    sxx = float(np.sum((x - x_bar) ** 2))
    slope = float(np.sum((x - x_bar) * (y - y_bar)) / sxx)
    intercept = y_bar - slope * x_bar
    y_fit = intercept + slope * x
    resid = y - y_fit
    sigma2 = float(np.sum(resid ** 2) / max(n - 2, 1))
    se = np.sqrt(np.maximum(0.0, sigma2 * (1.0 / n + ((x_grid - x_bar) ** 2) / sxx)))
    z = 1.96
    mid = intercept + slope * x_grid
    return mid, mid - z * se, mid + z * se


def style_panel_jk_axis(ax) -> None:
    ax.set_facecolor("#FCFCFC")
    ax.set_axisbelow(True)
    ax.grid(True, color="#DDDDDD", alpha=0.55, linewidth=0.6)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("#6C6C6C")
    ax.tick_params(labelsize=8, width=0.8, colors="#333333")


def plot_panel_j_reference(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    correct_mask: np.ndarray,
    title: str,
    show_ylabel: bool,
    show_xlabel: bool,
) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    correct_mask = np.asarray(correct_mask).astype(bool)
    xx = np.linspace(0.0, 1.0, 160)
    yy, lo, hi = linear_fit_ci_band(x, y, xx)
    rho = float(spearmanr(x, y).correlation) if x.size >= 3 else np.nan

    ax.fill_between(xx, lo, hi, color="#D7D7D7", alpha=0.38, linewidth=0, zorder=1)
    ax.plot(xx, yy, color="#2F6BFF", lw=2.2, zorder=4)
    ax.scatter(
        x[~correct_mask],
        y[~correct_mask],
        s=11,
        color="#555555",
        alpha=0.34,
        marker="o",
        linewidths=0,
        zorder=2,
    )
    ax.scatter(
        x[correct_mask],
        y[correct_mask],
        s=20,
        color="black",
        alpha=0.58,
        marker="^",
        linewidths=0,
        zorder=3,
    )

    y_pad = max(0.03, 0.06 * float(np.ptp(y)) if np.ptp(y) > 0 else 0.05)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(float(np.min(y)) - y_pad, float(np.max(y)) + y_pad)
    ax.set_title(title, fontsize=10, pad=6)
    if np.isfinite(rho):
        ax.text(
            0.04,
            0.95,
            f"rho={rho:.2f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            color="#4A4A4A",
        )
    if show_ylabel:
        ax.set_ylabel("Predicted Probability for the Positive Class", fontsize=9)
    if show_xlabel:
        ax.set_xlabel("Normalized Feature Values [0 1]", fontsize=9)
    style_panel_jk_axis(ax)
    ax.set_box_aspect(1.16)


def plot_panel_k_reference(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    correct_mask: np.ndarray,
    title: str,
    show_ylabel: bool,
    show_xlabel: bool,
) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    correct_mask = np.asarray(correct_mask).astype(bool)
    xx = np.linspace(0.0, 1.0, 160)
    yy, lo, hi = linear_fit_ci_band(x, y, xx)
    rho = float(spearmanr(x, y).correlation) if x.size >= 3 else np.nan

    ax.fill_between(xx, lo, hi, color="#D7D7D7", alpha=0.38, linewidth=0, zorder=1)
    ax.plot(xx, yy, color="#2F6BFF", lw=2.2, zorder=4)
    ax.scatter(
        x[~correct_mask],
        y[~correct_mask],
        s=11,
        color="#555555",
        alpha=0.34,
        marker="o",
        linewidths=0,
        zorder=2,
    )
    ax.scatter(
        x[correct_mask],
        y[correct_mask],
        s=20,
        color="black",
        alpha=0.58,
        marker="^",
        linewidths=0,
        zorder=3,
    )

    y_pad = max(0.025, 0.08 * float(np.ptp(y)) if np.ptp(y) > 0 else 0.05)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(float(np.min(y)) - y_pad, float(np.max(y)) + y_pad)
    ax.set_title(title, fontsize=10, pad=6)
    if np.isfinite(rho):
        ax.text(
            0.04,
            0.95,
            f"rho={rho:.2f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            color="#4A4A4A",
        )
    if show_ylabel:
        ax.set_ylabel("SHAP values", fontsize=9)
    if show_xlabel:
        ax.set_xlabel("normalized feature values [0 1]", fontsize=9)
    style_panel_jk_axis(ax)
    ax.set_box_aspect(1.16)


def load_layout_meta(path: str) -> Tuple[float, float, Dict[str, Dict[str, float]]]:
    with open(path, "r", encoding="utf-8") as f:
        obj = json.load(f)

    fig_obj = obj.get("figure", {})
    fig_w = float(fig_obj.get("w", 0))
    fig_h = float(fig_obj.get("h", 0))
    if fig_w <= 0 or fig_h <= 0:
        raise RuntimeError(f"Malformed layout metadata figure size: {path}")

    panels = obj.get("panels", [])
    panel_map: Dict[str, Dict[str, float]] = {}
    for row in panels:
        pid = str(row.get("panel_id", "")).strip()
        if not pid:
            continue
        panel_map[pid] = {
            "x": float(row["x"]),
            "y": float(row["y"]),
            "w": float(row["w"]),
            "h": float(row["h"]),
        }
    if not panel_map:
        raise RuntimeError(f"No panels found in layout metadata: {path}")
    return fig_w, fig_h, panel_map


def panel_rect_norm(panel: Dict[str, float], fig_w: float, fig_h: float) -> List[float]:
    left = panel["x"] / fig_w
    bottom = 1.0 - ((panel["y"] + panel["h"]) / fig_h)
    width = panel["w"] / fig_w
    height = panel["h"] / fig_h
    return [left, bottom, width, height]


def add_axes_from_panel(
    fig: plt.Figure,
    panel: Dict[str, float],
    fig_w: float,
    fig_h: float,
) -> plt.Axes:
    return fig.add_axes(panel_rect_norm(panel, fig_w, fig_h))


def add_gridspec_from_panel(
    fig: plt.Figure,
    panel: Dict[str, float],
    fig_w: float,
    fig_h: float,
    nrows: int,
    ncols: int,
    *,
    wspace: float = 0.2,
    hspace: float = 0.2,
    width_ratios: List[float] | None = None,
    height_ratios: List[float] | None = None,
):
    left, bottom, width, height = panel_rect_norm(panel, fig_w, fig_h)
    return fig.add_gridspec(
        nrows=nrows,
        ncols=ncols,
        left=left,
        right=left + width,
        bottom=bottom,
        top=bottom + height,
        wspace=wspace,
        hspace=hspace,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )


def require_layout_panels(
    panel_map: Dict[str, Dict[str, float]],
    names: List[str],
    layout_path: str,
) -> List[Dict[str, float]]:
    missing = [name for name in names if name not in panel_map]
    if missing:
        raise RuntimeError(f"Layout metadata missing panels {missing}: {layout_path}")
    return [panel_map[name] for name in names]


def fallback_coords_from_adata(adata: ad.AnnData, seed: int) -> Tuple[np.ndarray, str]:
    if "X_umap" in adata.obsm and adata.obsm["X_umap"].shape[1] >= 2:
        return np.asarray(adata.obsm["X_umap"][:, :2]), "X_umap"
    if "X_pca" in adata.obsm and adata.obsm["X_pca"].shape[1] >= 2:
        return np.asarray(adata.obsm["X_pca"][:, :2]), "X_pca"

    # Last-resort fallback: derive 2D coordinates from current expression matrix.
    from sklearn.decomposition import PCA

    mat = to_dense(adata.X)
    return PCA(n_components=2, random_state=seed).fit_transform(mat), "PCA(from X)"


def load_coords(adata: ad.AnnData, umap_csv: str, seed: int) -> Tuple[np.ndarray, str]:
    obs = np.asarray(adata.obs_names.astype(str))
    if umap_csv and os.path.exists(umap_csv):
        cdf = pd.read_csv(umap_csv)
        req = {"cell_id", "UMAP_1", "UMAP_2"}
        if req.issubset(set(cdf.columns)):
            cdf = cdf.drop_duplicates(subset=["cell_id"]).copy()
            cdf["cell_id"] = cdf["cell_id"].astype(str)
            cdf = cdf.set_index("cell_id")

            aligned = cdf.reindex(obs)[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
            ok = np.isfinite(aligned).all(axis=1)
            n_ok = int(np.sum(ok))
            if n_ok == adata.n_obs:
                return aligned, f"seurat_umap_csv({n_ok}/{adata.n_obs})"

            fb, fb_source = fallback_coords_from_adata(adata, seed)
            out = fb.copy()
            out[ok] = aligned[ok]
            n_miss = int(np.sum(~ok))
            src = f"seurat_umap_csv({n_ok}/{adata.n_obs}) + {fb_source}(missing={n_miss})"
            return out, src

    fb, fb_source = fallback_coords_from_adata(adata, seed)
    return fb, fb_source


def compute_umap_from_pca(adata: ad.AnnData, seed: int) -> Tuple[np.ndarray, str]:
    if "X_pca" not in adata.obsm or adata.obsm["X_pca"].shape[1] < 2:
        return fallback_coords_from_adata(adata, seed)

    import umap

    pcs = np.asarray(adata.obsm["X_pca"])
    n_pcs = min(30, pcs.shape[1])
    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.3,
        metric="euclidean",
        random_state=seed,
    )
    coords = reducer.fit_transform(pcs[:, :n_pcs])
    return np.asarray(coords), f"umap_from_X_pca(n_pcs={n_pcs})"


def load_or_compute_cached_umap_coords(
    adata: ad.AnnData, umap_csv: str, seed: int
) -> Tuple[np.ndarray, str]:
    obs = np.asarray(adata.obs_names.astype(str))
    if umap_csv and os.path.exists(umap_csv):
        cdf = pd.read_csv(umap_csv)
        req = {"cell_id", "UMAP_1", "UMAP_2"}
        if req.issubset(set(cdf.columns)):
            cdf = cdf.drop_duplicates(subset=["cell_id"]).copy()
            cdf["cell_id"] = cdf["cell_id"].astype(str)
            cdf = cdf.set_index("cell_id")
            aligned = cdf.reindex(obs)[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
            if np.isfinite(aligned).all():
                return aligned, f"cached_umap_csv({adata.n_obs}/{adata.n_obs})"

    coords, source = compute_umap_from_pca(adata, seed)
    if umap_csv:
        out_df = pd.DataFrame(
            {
                "cell_id": obs,
                "UMAP_1": coords[:, 0],
                "UMAP_2": coords[:, 1],
            }
        )
        out_df.to_csv(umap_csv, index=False)
    return coords, source


def r_quote(path: str) -> str:
    return path.replace("\\", "/").replace("'", "\\'")


def compute_malignant_umap_from_harmony(
    cell_ids: np.ndarray,
    seurat_rds: str,
    out_csv: str,
    seed: int,
) -> Tuple[str, str]:
    if not seurat_rds or not os.path.exists(seurat_rds):
        raise RuntimeError(f"Missing Seurat RDS for malignant Harmony UMAP: {seurat_rds}")

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="fig6_malignant_umap_") as tmpdir:
        cell_csv = os.path.join(tmpdir, "cell_ids.csv")
        pd.DataFrame({"cell_id": np.asarray(cell_ids, dtype=str)}).to_csv(cell_csv, index=False)

        r_script = f"""
suppressPackageStartupMessages({{library(Seurat); library(data.table); library(uwot)}})
seu <- readRDS('{r_quote(seurat_rds)}')
cells <- fread('{r_quote(cell_csv)}', data.table = FALSE)$cell_id
emb <- Embeddings(seu, 'harmony')
keep <- intersect(as.character(cells), rownames(emb))
if (length(keep) < 100) stop('Too few malignant cells found in Harmony embeddings: ', length(keep))
hm <- emb[match(keep, rownames(emb)), , drop = FALSE]
set.seed({int(seed)})
um <- uwot::umap(
  hm,
  n_neighbors = 30,
  min_dist = 0.3,
  metric = 'cosine',
  n_threads = {{
    ncpu <- suppressWarnings(as.integer(parallel::detectCores()))
    if (is.na(ncpu) || ncpu < 2) 1L else (ncpu - 1L)
  }},
  ret_model = FALSE,
  verbose = FALSE
)
out <- data.frame(cell_id = keep, UMAP_1 = um[, 1], UMAP_2 = um[, 2], stringsAsFactors = FALSE)
write.csv(out, '{r_quote(out_csv)}', row.names = FALSE)
cat('saved', '{r_quote(out_csv)}', '\\n')
"""
        proc = subprocess.run(
            ["Rscript", "-"],
            input=r_script,
            text=True,
            capture_output=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "Failed to compute malignant Harmony UMAP via Rscript.\n"
                f"stdout:\n{proc.stdout}\n\nstderr:\n{proc.stderr}"
            )
    return out_csv, "harmony_umap_rscript"


def load_or_compute_malignant_harmony_coords(
    cell_ids: np.ndarray,
    umap_csv: str,
    seurat_rds: str,
    seed: int,
    fallback_coords: np.ndarray | None = None,
    fallback_source: str = "",
) -> Tuple[np.ndarray, str]:
    obs = pd.Index(np.asarray(cell_ids, dtype=str), name="cell_id")

    def load_cached(path: str) -> Tuple[np.ndarray, str]:
        cdf = pd.read_csv(path)
        req = {"cell_id", "UMAP_1", "UMAP_2"}
        if not req.issubset(set(cdf.columns)):
            raise RuntimeError(f"Malformed malignant UMAP CSV: {path}")
        cdf = cdf.drop_duplicates(subset=["cell_id"]).copy()
        cdf["cell_id"] = cdf["cell_id"].astype(str)
        cdf = cdf.set_index("cell_id")
        aligned = cdf.reindex(obs)[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
        ok = np.isfinite(aligned).all(axis=1)
        if not np.all(ok):
            missing = int(np.sum(~ok))
            if fallback_coords is None or fallback_coords.shape != aligned.shape:
                raise RuntimeError(f"Malignant UMAP CSV missing {missing} cells: {path}")
            aligned[~ok] = fallback_coords[~ok]
            src = f"cached_harmony_umap_csv({int(np.sum(ok))}/{len(obs)}) + {fallback_source}(missing={missing})"
            return aligned, src
        return aligned, f"cached_harmony_umap_csv({len(obs)}/{len(obs)})"

    if umap_csv and os.path.exists(umap_csv):
        try:
            return load_cached(umap_csv)
        except Exception:
            pass

    if not umap_csv:
        raise RuntimeError("Panel B requires a writable malignant UMAP CSV path.")
    out_csv, source = compute_malignant_umap_from_harmony(cell_ids, seurat_rds, umap_csv, seed)
    coords, src = load_cached(out_csv)
    if src.startswith("cached_harmony_umap_csv("):
        src = source if f"({len(obs)}/{len(obs)})" in src else src.replace("cached_harmony_umap_csv", source)
    return coords, src


def maybe_load_cached_malignant_harmony_coords(
    cell_ids: np.ndarray,
    umap_csv: str,
) -> Tuple[np.ndarray | None, str]:
    obs = pd.Index(np.asarray(cell_ids, dtype=str), name="cell_id")
    if not umap_csv or not os.path.exists(umap_csv):
        return None, ""
    cdf = pd.read_csv(umap_csv)
    req = {"cell_id", "UMAP_1", "UMAP_2"}
    if not req.issubset(set(cdf.columns)):
        return None, ""
    cdf = cdf.drop_duplicates(subset=["cell_id"]).copy()
    cdf["cell_id"] = cdf["cell_id"].astype(str)
    cdf = cdf.set_index("cell_id")
    aligned = cdf.reindex(obs)[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)
    if np.isfinite(aligned).all():
        return aligned, f"cached_harmony_umap_csv({len(obs)}/{len(obs)})"
    return None, ""


def polygon_area(poly: np.ndarray) -> float:
    if poly.shape[0] < 3:
        return 0.0
    x = poly[:, 0]
    y = poly[:, 1]
    return 0.5 * abs(np.dot(x[:-1], y[1:]) - np.dot(y[:-1], x[1:]))


def convex_hull_polygon(points: np.ndarray) -> np.ndarray | None:
    if points.shape[0] < 3:
        return None
    try:
        hull = ConvexHull(points)
    except Exception:
        return None
    verts = hull.vertices.tolist()
    poly = points[verts]
    if poly.shape[0] < 3:
        return None
    return np.vstack([poly, poly[0]])


def triangle_circumradius(tri_pts: np.ndarray) -> np.ndarray:
    a = np.linalg.norm(tri_pts[:, 1] - tri_pts[:, 0], axis=1)
    b = np.linalg.norm(tri_pts[:, 2] - tri_pts[:, 1], axis=1)
    c = np.linalg.norm(tri_pts[:, 0] - tri_pts[:, 2], axis=1)
    area2 = np.abs(
        (tri_pts[:, 1, 0] - tri_pts[:, 0, 0]) * (tri_pts[:, 2, 1] - tri_pts[:, 0, 1])
        - (tri_pts[:, 1, 1] - tri_pts[:, 0, 1]) * (tri_pts[:, 2, 0] - tri_pts[:, 0, 0])
    )
    area = area2 / 2.0
    out = np.full(tri_pts.shape[0], np.inf, dtype=float)
    ok = area > 1e-12
    out[ok] = (a[ok] * b[ok] * c[ok]) / (4.0 * area[ok])
    return out


def connected_components(nodes: List[int], adjacency: Dict[int, set[int]]) -> List[List[int]]:
    out: List[List[int]] = []
    seen: set[int] = set()
    for node in nodes:
        if node in seen:
            continue
        stack = [node]
        comp: List[int] = []
        seen.add(node)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in adjacency[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        out.append(comp)
    return out


def trace_loop(component_nodes: List[int], adjacency: Dict[int, set[int]], points: np.ndarray) -> np.ndarray | None:
    degrees = [len(adjacency[n]) for n in component_nodes]
    if not degrees or any(d != 2 for d in degrees):
        return None

    start = min(component_nodes, key=lambda idx: (points[idx, 0], points[idx, 1]))
    centroid = points[component_nodes].mean(axis=0)
    nbrs = list(adjacency[start])
    if len(nbrs) != 2:
        return None
    ang = [np.arctan2(points[n, 1] - centroid[1], points[n, 0] - centroid[0]) for n in nbrs]
    curr = nbrs[int(np.argmin(ang))]
    prev = start
    order = [start]

    for _ in range(len(component_nodes) + 2):
        order.append(curr)
        nxts = [n for n in adjacency[curr] if n != prev]
        if not nxts:
            return None
        nxt = nxts[0]
        prev, curr = curr, nxt
        if curr == start:
            break
    if order[-1] != start:
        order.append(start)
    poly = points[order]
    return poly if poly.shape[0] >= 4 else None


def alpha_shape_polygons(points: np.ndarray, seed: int) -> List[np.ndarray]:
    if points.shape[0] < 4:
        poly = convex_hull_polygon(points)
        return [poly] if poly is not None else []

    pts = np.unique(np.asarray(points, dtype=float), axis=0)
    if pts.shape[0] < 4:
        poly = convex_hull_polygon(pts)
        return [poly] if poly is not None else []

    max_points = 12000
    if pts.shape[0] > max_points:
        rng = np.random.default_rng(seed)
        keep = rng.choice(pts.shape[0], size=max_points, replace=False)
        pts = pts[keep]

    try:
        tree = cKDTree(pts)
        k = min(9, max(3, pts.shape[0] - 1))
        dists, _ = tree.query(pts, k=k + 1)
        local_scale = float(np.median(dists[:, -1]))
        tri = Delaunay(pts)
    except Exception:
        poly = convex_hull_polygon(pts)
        return [poly] if poly is not None else []

    tri_pts = pts[tri.simplices]
    radii = triangle_circumradius(tri_pts)
    radius_cut = max(local_scale * 2.2, 1e-6)
    keep = tri.simplices[np.isfinite(radii) & (radii <= radius_cut)]
    if keep.shape[0] == 0:
        poly = convex_hull_polygon(pts)
        return [poly] if poly is not None else []

    edge_count: Dict[Tuple[int, int], int] = {}
    for a, b, c in keep:
        for i, j in ((a, b), (b, c), (c, a)):
            key = (i, j) if i < j else (j, i)
            edge_count[key] = edge_count.get(key, 0) + 1

    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
    if len(boundary_edges) < 3:
        poly = convex_hull_polygon(pts)
        return [poly] if poly is not None else []

    adjacency: Dict[int, set[int]] = {}
    for i, j in boundary_edges:
        adjacency.setdefault(i, set()).add(j)
        adjacency.setdefault(j, set()).add(i)

    polys: List[np.ndarray] = []
    components = connected_components(list(adjacency.keys()), adjacency)
    full_area = float((pts[:, 0].max() - pts[:, 0].min()) * (pts[:, 1].max() - pts[:, 1].min()))
    min_area = max(full_area * 0.0015, 1e-4)
    for comp in components:
        if len(comp) < 3:
            continue
        poly = trace_loop(comp, adjacency, pts)
        if poly is None:
            poly = convex_hull_polygon(pts[np.asarray(comp)])
        if poly is None or poly.shape[0] < 4:
            continue
        if polygon_area(poly) >= min_area:
            polys.append(poly)

    if polys:
        polys.sort(key=polygon_area, reverse=True)
        return polys

    poly = convex_hull_polygon(pts)
    return [poly] if poly is not None else []


def compute_frame_and_outline(
    coords: np.ndarray,
    seed: int,
) -> Tuple[Tuple[float, float], Tuple[float, float], List[np.ndarray] | None]:
    valid = np.isfinite(coords).all(axis=1)
    c = coords[valid]
    if c.shape[0] < 3:
        return (-1.0, 1.0), (-1.0, 1.0), None

    outline = convex_hull_polygon(c)
    x_min, x_max = float(np.min(c[:, 0])), float(np.max(c[:, 0]))
    y_min, y_max = float(np.min(c[:, 1])), float(np.max(c[:, 1]))
    # Keep a visible black margin between the hull and panel edge, like the paper.
    x_pad = max((x_max - x_min) * 0.055, 1e-6)
    y_pad = max((y_max - y_min) * 0.055, 1e-6)
    return (x_min - x_pad, x_max + x_pad), (y_min - y_pad, y_max + y_pad), [outline] if outline is not None else None


def save_abc_preview(
    out_dir: str,
    genes: List[str],
    expr: np.ndarray,
    expr_malignant: np.ndarray,
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    xlim_a: Tuple[float, float],
    ylim_a: Tuple[float, float],
    xlim_b: Tuple[float, float],
    ylim_b: Tuple[float, float],
    outline_a: List[np.ndarray] | None,
    outline_b: List[np.ndarray] | None,
    y_c: np.ndarray,
    dpi: int,
) -> Tuple[str, str]:
    n_genes = len(genes)
    fig = plt.figure(figsize=(24.0, 9.2), facecolor="#EBEBEB")
    outer = GridSpec(
        nrows=2,
        ncols=3,
        figure=fig,
        height_ratios=[1.0, 1.0],
        width_ratios=[2.0, 1.45, 1.20],
        hspace=0.26,
        wspace=0.20,
    )
    fig.subplots_adjust(left=0.045, right=0.975, top=0.93, bottom=0.10)

    gs_a = outer[0, 0:2].subgridspec(1, n_genes, wspace=0.20)
    axes_a = [fig.add_subplot(gs_a[0, i]) for i in range(n_genes)]
    for i, g in enumerate(genes):
        weighted_density(
            axes_a[i],
            coords_a,
            expr[:, i],
            g,
            xlim=xlim_a,
            ylim=ylim_a,
            outline_xy=outline_a,
        )
    add_panel_letter(axes_a[0], "A")

    ax_c = fig.add_subplot(outer[0, 2])
    colors = sns.color_palette("tab10", n_colors=n_genes)
    for i, g in enumerate(genes):
        score = expr[:, i]
        fpr, tpr, _ = roc_curve(y_c, score)
        auc_i = auc(fpr, tpr)
        ax_c.plot(fpr, tpr, lw=2.5, color=colors[i], label=f"{g}, AUC={auc_i:.3f}")
    ax_c.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax_c.set_xlabel("1 - Specificity", fontsize=18)
    ax_c.set_ylabel("Sensitivity", fontsize=18)
    ax_c.tick_params(labelsize=13)
    ax_c.legend(loc="lower right", fontsize=12, frameon=False)
    add_panel_letter(ax_c, "C")

    gs_b = outer[1, 0:2].subgridspec(1, n_genes, wspace=0.20)
    axes_b = [fig.add_subplot(gs_b[0, i]) for i in range(n_genes)]
    for i, g in enumerate(genes):
        weighted_density(
            axes_b[i],
            coords_b,
            expr_malignant[:, i],
            g,
            xlim=xlim_b,
            ylim=ylim_b,
            outline_xy=outline_b,
        )
    add_panel_letter(axes_b[0], "B")

    ax_blank = fig.add_subplot(outer[1, 2])
    ax_blank.axis("off")

    ax_c.set_facecolor("#EBEBEB")
    ax_blank.set_facecolor("#EBEBEB")

    out_png = os.path.join(out_dir, "Fig6_preview_ABC.png")
    out_pdf = os.path.join(out_dir, "Fig6_preview_ABC.pdf")
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return out_png, out_pdf


def save_def_preview(
    out_dir: str,
    fit_df: pd.DataFrame,
    fits: List[FitResult],
    best_fit: FitResult,
    y_train: np.ndarray,
    y_test: np.ndarray,
    dpi: int,
) -> Tuple[str, str]:
    fig = plt.figure(figsize=(16.0, 5.8), facecolor="#EBEBEB")
    outer = GridSpec(
        nrows=1,
        ncols=3,
        figure=fig,
        width_ratios=[1.50, 1.35, 1.18],
        wspace=0.32,
    )
    fig.subplots_adjust(left=0.055, right=0.985, top=0.92, bottom=0.20)

    model_order = figure6_model_order([f.name for f in fits])
    model_palette = figure6_model_palette(model_order)

    ax_d = fig.add_subplot(outer[0, 0])
    sns.boxplot(
        data=fit_df,
        x="model",
        y="auc",
        order=model_order,
        hue="model",
        palette=model_palette,
        dodge=False,
        ax=ax_d,
        fliersize=1.8,
        linewidth=1.0,
        saturation=0.95,
    )
    if ax_d.legend_ is not None:
        ax_d.legend_.remove()
    sns.stripplot(
        data=fit_df,
        x="model",
        y="auc",
        order=model_order,
        ax=ax_d,
        color="#555555",
        size=2.0,
        alpha=0.55,
        jitter=0.15,
    )
    ax_d.set_title("Machine Learning model benchmark", fontsize=13, fontweight="bold")
    ax_d.set_xlabel("Learn Model", fontsize=9)
    ax_d.set_ylabel("Mean AUC", fontsize=9)
    ax_d.tick_params(axis="x", rotation=42, labelsize=8)
    ax_d.tick_params(axis="y", labelsize=8)
    ax_d.grid(axis="y", alpha=0.25, color="#CFCFCF")
    add_panel_letter(ax_d, "D")

    gs_e = outer[0, 1].subgridspec(1, 2, wspace=0.28)
    ax_e1 = fig.add_subplot(gs_e[0, 0])
    ax_e2 = fig.add_subplot(gs_e[0, 1])
    for name in model_order:
        fit = next(f for f in fits if f.name == name)
        c = model_palette[name]
        prec, rec, _ = precision_recall_curve(y_test, fit.test_proba)
        fpr, tpr, _ = roc_curve(y_test, fit.test_proba)
        ax_e1.plot(rec, prec, lw=1.8, color=c, label=name)
        ax_e2.plot(fpr, tpr, lw=1.8, color=c, label=name)
    ax_e1.set_xlim(0, 1)
    ax_e1.set_ylim(0.5, 1.01)
    ax_e1.set_xlabel("Recall", fontsize=9)
    ax_e1.set_ylabel("Precision", fontsize=9)
    ax_e1.tick_params(labelsize=8)
    ax_e2.set_xlim(0, 1)
    ax_e2.set_ylim(0, 1.01)
    ax_e2.set_xlabel("1 - Specificity", fontsize=9)
    ax_e2.set_ylabel("Sensitivity", fontsize=9)
    ax_e2.tick_params(labelsize=8)
    ax_e2.plot([0, 1], [0, 1], "--", color="lightgray", lw=1)
    ax_e2.legend(
        title="Learner",
        fontsize=7,
        title_fontsize=8,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )
    add_panel_letter(ax_e1, "E")

    gs_f = outer[0, 2].subgridspec(2, 2, wspace=0.36, hspace=0.42)
    ax_f11 = fig.add_subplot(gs_f[0, 0])
    ax_f12 = fig.add_subplot(gs_f[0, 1])
    ax_f21 = fig.add_subplot(gs_f[1, 0])
    ax_f22 = fig.add_subplot(gs_f[1, 1])

    fpr_tr, tpr_tr, _ = roc_curve(y_train, best_fit.train_proba)
    fpr_te, tpr_te, _ = roc_curve(y_test, best_fit.test_proba)
    auc_tr = auc(fpr_tr, tpr_tr)
    auc_te = auc(fpr_te, tpr_te)
    prec_tr, rec_tr, _ = precision_recall_curve(y_train, best_fit.train_proba)
    prec_te, rec_te, _ = precision_recall_curve(y_test, best_fit.test_proba)
    pr_tr = auc(rec_tr, prec_tr)
    pr_te = auc(rec_te, prec_te)
    line_c = "#7FD0D5"
    base_train = float(np.mean(y_train))
    base_test = float(np.mean(y_test))

    ax_f11.plot(fpr_tr, tpr_tr, color=line_c, lw=2)
    ax_f11.plot([0, 1], [0, 1], "--", color="lightgray", lw=1)
    ax_f11.set_title(f"development set N={len(y_train)}, AUC = {auc_tr:.2f}", fontsize=6.3, pad=2)
    ax_f11.set_xlabel("False Positive Rate (FPR)", fontsize=7)
    ax_f11.set_ylabel("True Positive Rate (TPR)", fontsize=7)
    ax_f12.plot(rec_tr, prec_tr, color=line_c, lw=2)
    ax_f12.axhline(base_train, color="lightgray", lw=0.9, ls="--")
    ax_f12.set_title(f"development set N={len(y_train)}, PRAUC = {pr_tr:.2f}", fontsize=6.3, pad=2)
    ax_f12.set_xlabel("Recall", fontsize=7)
    ax_f12.set_ylabel("Precision", fontsize=7)
    ax_f21.plot(fpr_te, tpr_te, color=line_c, lw=2)
    ax_f21.plot([0, 1], [0, 1], "--", color="lightgray", lw=1)
    ax_f21.set_title(f"test set N={len(y_test)}, AUC = {auc_te:.2f}", fontsize=6.3, pad=2)
    ax_f21.set_xlabel("False Positive Rate (FPR)", fontsize=7)
    ax_f21.set_ylabel("True Positive Rate (TPR)", fontsize=7)
    ax_f22.plot(rec_te, prec_te, color=line_c, lw=2)
    ax_f22.axhline(base_test, color="lightgray", lw=0.9, ls="--")
    ax_f22.set_title(f"test set N={len(y_test)}, PRAUC = {pr_te:.2f}", fontsize=6.3, pad=2)
    ax_f22.set_xlabel("Recall", fontsize=7)
    ax_f22.set_ylabel("Precision", fontsize=7)
    for ax in [ax_f11, ax_f12, ax_f21, ax_f22]:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.tick_params(labelsize=7)
    add_panel_letter(ax_f11, "F")

    for ax in [ax_d, ax_e1, ax_e2, ax_f11, ax_f12, ax_f21, ax_f22]:
        ax.set_facecolor("#EBEBEB")

    out_png = os.path.join(out_dir, "Fig6_preview_DEF.png")
    out_pdf = os.path.join(out_dir, "Fig6_preview_DEF.pdf")
    fig.savefig(out_png, dpi=dpi, facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, facecolor=fig.get_facecolor())
    plt.close(fig)
    return out_png, out_pdf


def main() -> int:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.environ.get("MPLCONFIGDIR", "/tmp/mpl"), exist_ok=True)

    adata = ad.read_h5ad(args.adata_path)
    meta = read_metadata_csv(args.meta_path)
    panel_b_meta_path = args.meta_path_panel_b or infer_panel_b_meta_path(args.run_dir) or args.meta_path
    if not os.path.exists(panel_b_meta_path):
        raise RuntimeError(f"Missing panel B metadata CSV: {panel_b_meta_path}")
    meta_panel_b = read_metadata_csv(panel_b_meta_path)
    classifier_meta_path = args.meta_path_classifier or os.path.join(
        args.run_dir, "tables", "metadata_strict_latest13k.csv"
    )
    if not os.path.exists(classifier_meta_path):
        raise RuntimeError(f"Missing classifier metadata CSV: {classifier_meta_path}")
    meta_classifier = read_metadata_csv(classifier_meta_path)

    for col in ["stem_group_strict"]:
        if col not in meta.columns:
            raise RuntimeError(f"Missing metadata column: {col}")
    if "stem_group_strict" not in meta_classifier.columns:
        raise RuntimeError("Missing stem_group_strict in classifier metadata CSV.")
    if "malignant_status" not in adata.obs.columns:
        raise RuntimeError("Missing adata.obs column: malignant_status")

    adata_cells = np.asarray(adata.obs_names.astype(str))
    coords_a, coord_source_a = load_coords(adata, args.umap_csv, args.seed)

    gene_file = args.gene_file or os.path.join(args.run_dir, "tables", "final_consensus_genes.csv")
    genes = load_genes(gene_file, adata)
    n_genes = len(genes)
    if n_genes < 2:
        raise RuntimeError("Too few genes for Figure 6.")

    if adata.raw is not None and all(g in set(map(str, adata.raw.var_names)) for g in genes):
        expr = to_dense(adata.raw[:, genes].X)
    else:
        missing = [g for g in genes if g not in set(map(str, adata.var_names))]
        if missing:
            raise RuntimeError(f"Genes missing from adata: {missing}")
        expr = to_dense(adata[:, genes].X)

    # Panel A uses the full atlas. Panel B and downstream panels use their own aligned metadata.
    meta_index = set(meta.index.astype(str))
    common = [c for c in adata_cells if c in meta_index]
    if not common:
        raise RuntimeError("No common cells between adata and metadata_strict.")
    adata_common = adata[np.asarray(common), :].copy()
    meta = meta.loc[common].copy()

    panel_b_index = set(meta_panel_b.index.astype(str))
    common_b = [c for c in adata_cells if c in panel_b_index]
    if not common_b:
        raise RuntimeError("No common cells between adata and panel B metadata.")
    adata_panel_b = adata[np.asarray(common_b), :].copy()
    meta_panel_b = meta_panel_b.loc[common_b].copy()

    if "malignant_status" in meta_panel_b.columns:
        malignant_meta_col_b = "malignant_status"
    elif "malignant_status_recalc" in meta_panel_b.columns:
        malignant_meta_col_b = "malignant_status_recalc"
    else:
        raise RuntimeError("Missing malignant status column in panel B metadata.")

    malignant_mask_b = meta_panel_b[malignant_meta_col_b].astype(str).values == "malignant"
    if int(np.sum(malignant_mask_b)) < 50:
        raise RuntimeError("Too few malignant cells after panel B filtering.")

    adata_malignant = adata_panel_b[malignant_mask_b, :].copy()
    meta_malignant = meta_panel_b.loc[malignant_mask_b].copy()

    if adata_malignant.raw is not None and all(
        g in set(map(str, adata_malignant.raw.var_names)) for g in genes
    ):
        expr_malignant = to_dense(adata_malignant.raw[:, genes].X)
    else:
        expr_malignant = to_dense(adata_malignant[:, genes].X)

    umap_csv_b = panel_b_umap_cache_path(panel_b_meta_path, args.umap_csv_malignant)
    coords_b_cached, coord_source_b_cached = maybe_load_cached_malignant_harmony_coords(
        adata_malignant.obs_names.astype(str).to_numpy(), umap_csv_b
    )
    if coords_b_cached is not None:
        coords_b, coord_source_b = coords_b_cached, coord_source_b_cached
    else:
        coords_b_fallback, coord_source_b_fallback = compute_umap_from_pca(adata_malignant, args.seed)
        coords_b, coord_source_b = load_or_compute_malignant_harmony_coords(
            adata_malignant.obs_names.astype(str).to_numpy(),
            umap_csv_b,
            args.seurat_rds,
            args.seed,
            fallback_coords=coords_b_fallback,
            fallback_source=coord_source_b_fallback,
        )

    xlim_a, ylim_a, outline_a = compute_frame_and_outline(coords_a, seed=args.seed)
    xlim_b, ylim_b, outline_b = compute_frame_and_outline(coords_b, seed=args.seed)

    # Panel C follows the paper: single-gene ROC for malignant vs non-malignant cells.
    # copy: label panel C with the same 13,664 inferCNV malignant cells used in panel B
    y_c = adata.obs_names.astype(str).isin(set(adata_malignant.obs_names.astype(str))).astype(int)
    if len(np.unique(y_c)) < 2:
        raise RuntimeError("Panel C requires both malignant and non-malignant cells.")

    if args.preview_only_abc:
        out_png, out_pdf = save_abc_preview(
            out_dir=args.out_dir,
            genes=genes,
            expr=expr,
            expr_malignant=expr_malignant,
            coords_a=coords_a,
            coords_b=coords_b,
            xlim_a=xlim_a,
            ylim_a=ylim_a,
            xlim_b=xlim_b,
            ylim_b=ylim_b,
            outline_a=outline_a,
            outline_b=outline_b,
            y_c=y_c,
            dpi=args.dpi,
        )
        print(f"[OK] Saved: {out_png}")
        print(f"[OK] Saved: {out_pdf}")
        print(f"[INFO] Preview only: A/B/C")
        print(f"[INFO] Panel A cells(all atlas): {adata.n_obs}")
        print(f"[INFO] Panel B cells(malignant subset): {adata_malignant.n_obs}")
        print(f"[INFO] Coordinate source A: {coord_source_a}")
        print(f"[INFO] Coordinate source B: {coord_source_b}")
        print(f"[INFO] Panel B metadata: {panel_b_meta_path}")
        print(f"[INFO] Panel C label counts: malignant={int(np.sum(y_c == 1))}, non_malignant={int(np.sum(y_c == 0))}")
        print(f"[INFO] Gene file: {gene_file}")
        print(f"[INFO] Genes used ({n_genes}): {genes}")
        return 0

    # Panels D-K follow the run-specific strict cohort used during upstream gene screening.
    cls_common = [c for c in adata_cells if c in set(meta_classifier.index.astype(str))]
    if not cls_common:
        raise RuntimeError("No common cells between adata and classifier metadata CSV.")
    adata_classifier = adata[np.asarray(cls_common), :].copy()
    meta_classifier = meta_classifier.loc[cls_common].copy()
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
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, random_state=args.seed, stratify=y
    )
    analysis_cache_path = args.analysis_cache or default_analysis_cache_path(
        args.run_dir, classifier_meta_path, genes, args.seed
    )
    cache_loaded = False
    fits: List[FitResult] = []

    if os.path.exists(analysis_cache_path) and not args.refresh_analysis_cache:
        try:
            cache = load_analysis_cache(analysis_cache_path)
            if (
                cache.get("cache_version") == 1
                and cache.get("genes") == genes
                and int(cache.get("seed", -1)) == int(args.seed)
            ):
                X_test = np.asarray(cache["X_test"], dtype=float)
                y_train = np.asarray(cache["y_train"], dtype=int)
                y_test = np.asarray(cache["y_test"], dtype=int)
                shap_vals = np.asarray(cache["shap_vals"], dtype=float)
                for rec in cache["fit_records"]:
                    fits.append(
                        FitResult(
                            name=rec["name"],
                            model=None,
                            cv_auc=[float(v) for v in rec["cv_auc"]],
                            train_proba=np.asarray(rec["train_proba"], dtype=float),
                            test_proba=np.asarray(rec["test_proba"], dtype=float),
                        )
                    )
                best_fit = next(f for f in fits if f.name == str(cache["best_name"]))
                cache_loaded = True
        except Exception:
            cache_loaded = False

    if not cache_loaded:
        models = model_zoo(args.seed)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=args.seed)
        for name, model in models.items():
            fold_auc: List[float] = []
            for tr, va in cv.split(X_train, y_train):
                model.fit(X_train[tr], y_train[tr])
                p = model.predict_proba(X_train[va])[:, 1]
                fold_auc.append(float(roc_auc_score(y_train[va], p)))
            model.fit(X_train, y_train)
            train_proba = model.predict_proba(X_train)[:, 1]
            test_proba = model.predict_proba(X_test)[:, 1]
            fits.append(
                FitResult(
                    name=name,
                    model=model,
                    cv_auc=fold_auc,
                    train_proba=train_proba,
                    test_proba=test_proba,
                )
            )

        best_fit = sorted(fits, key=lambda x: np.mean(x.cv_auc), reverse=True)[0]

        bg_n = min(120, X_train.shape[0])
        rng = np.random.default_rng(args.seed)
        bg_idx = rng.choice(X_train.shape[0], size=bg_n, replace=False)
        background = X_train[bg_idx]

        def best_predict_fn(data):
            arr = np.asarray(data, dtype=float)
            return best_fit.model.predict_proba(arr)[:, 1]

        shap_exp = shap.Explainer(best_predict_fn, background)(X_test)
        shap_vals = np.asarray(shap_exp.values)
        if shap_vals.ndim == 3 and shap_vals.shape[-1] == 1:
            shap_vals = shap_vals[:, :, 0]
        if shap_vals.ndim != 2:
            raise RuntimeError(f"Unexpected SHAP shape from best model: {shap_vals.shape}")

        save_analysis_cache(
            analysis_cache_path,
            {
                "cache_version": 1,
                "genes": genes,
                "seed": int(args.seed),
                "X_test": np.asarray(X_test, dtype=float),
                "y_train": np.asarray(y_train, dtype=int),
                "y_test": np.asarray(y_test, dtype=int),
                "shap_vals": np.asarray(shap_vals, dtype=float),
                "best_name": best_fit.name,
                "fit_records": [
                    {
                        "name": f.name,
                        "cv_auc": [float(v) for v in f.cv_auc],
                        "train_proba": np.asarray(f.train_proba, dtype=float),
                        "test_proba": np.asarray(f.test_proba, dtype=float),
                    }
                    for f in fits
                ],
            },
        )

    fit_df = pd.DataFrame(
        [(f.name, i + 1, a) for f in fits for i, a in enumerate(f.cv_auc)],
        columns=["model", "fold", "auc"],
    )
    model_order = figure6_model_order([f.name for f in fits])
    model_palette = figure6_model_palette(model_order)
    y_pred_train = (best_fit.train_proba >= 0.5).astype(int)
    y_pred_test = (best_fit.test_proba >= 0.5).astype(int)

    scaler_mm = MinMaxScaler()
    X_test_mm = scaler_mm.fit_transform(X_test)
    correct_test = (y_pred_test == y_test)

    if args.preview_only_def:
        out_png, out_pdf = save_def_preview(
            out_dir=args.out_dir,
            fit_df=fit_df,
            fits=fits,
            best_fit=best_fit,
            y_train=y_train,
            y_test=y_test,
            dpi=args.dpi,
        )
        print(f"[OK] Saved: {out_png}")
        print(f"[OK] Saved: {out_pdf}")
        print(f"[INFO] Preview only: D/E/F")
        print(f"[INFO] Analysis cache: {analysis_cache_path}")
        print(f"[INFO] Analysis cache loaded: {cache_loaded}")
        print(f"[INFO] Gene file: {gene_file}")
        print(f"[INFO] Genes used ({n_genes}): {genes}")
        print(f"[INFO] Train N={len(y_train)}, Test N={len(y_test)}, Best model={best_fit.name}")
        return 0

    layout_meta_path = args.layout_meta.strip()
    if layout_meta_path:
        if not os.path.exists(layout_meta_path):
            raise RuntimeError(f"Missing layout metadata JSON: {layout_meta_path}")
        layout_fig_w, layout_fig_h, panel_layout = load_layout_meta(layout_meta_path)
        if abs(args.fig_height - 17.5) < 1e-9:
            fig_height = args.fig_width * (layout_fig_h / layout_fig_w)
        else:
            fig_height = args.fig_height
        fig = plt.figure(figsize=(args.fig_width, fig_height), facecolor="#EBEBEB")

        panels_a = require_layout_panels(
            panel_layout,
            [f"A{i + 1}" for i in range(n_genes)],
            layout_meta_path,
        )
        panels_b = require_layout_panels(
            panel_layout,
            [f"B{i + 1}" for i in range(n_genes)],
            layout_meta_path,
        )
        panel_c = require_layout_panels(panel_layout, ["C"], layout_meta_path)[0]
        panel_d = require_layout_panels(panel_layout, ["D"], layout_meta_path)[0]
        panel_e = require_layout_panels(panel_layout, ["E"], layout_meta_path)[0]
        panel_f = require_layout_panels(panel_layout, ["F"], layout_meta_path)[0]
        panel_g = require_layout_panels(panel_layout, ["G"], layout_meta_path)[0]
        panel_h = require_layout_panels(panel_layout, ["H"], layout_meta_path)[0]
        panel_i = require_layout_panels(panel_layout, ["I"], layout_meta_path)[0]
        panel_j = require_layout_panels(panel_layout, ["J"], layout_meta_path)[0]
        panel_k = require_layout_panels(panel_layout, ["K"], layout_meta_path)[0]

        axes_a = [add_axes_from_panel(fig, panel, layout_fig_w, layout_fig_h) for panel in panels_a]
        ax_c = add_axes_from_panel(fig, panel_c, layout_fig_w, layout_fig_h)
        axes_b = [add_axes_from_panel(fig, panel, layout_fig_w, layout_fig_h) for panel in panels_b]
        ax_d = add_axes_from_panel(fig, panel_d, layout_fig_w, layout_fig_h)
        gs_e = add_gridspec_from_panel(
            fig, panel_e, layout_fig_w, layout_fig_h, 1, 2, wspace=0.28, hspace=0.0
        )
        ax_e1 = fig.add_subplot(gs_e[0, 0])
        ax_e2 = fig.add_subplot(gs_e[0, 1])
        gs_f = add_gridspec_from_panel(
            fig, panel_f, layout_fig_w, layout_fig_h, 2, 2, wspace=0.36, hspace=0.42
        )
        ax_f11 = fig.add_subplot(gs_f[0, 0])
        ax_f12 = fig.add_subplot(gs_f[0, 1])
        ax_f21 = fig.add_subplot(gs_f[1, 0])
        ax_f22 = fig.add_subplot(gs_f[1, 1])
        gs_g = add_gridspec_from_panel(
            fig, panel_g, layout_fig_w, layout_fig_h, 1, 2, wspace=0.25, hspace=0.0
        )
        ax_g1 = fig.add_subplot(gs_g[0, 0])
        ax_g2 = fig.add_subplot(gs_g[0, 1])
        ax_h = add_axes_from_panel(fig, panel_h, layout_fig_w, layout_fig_h)
        ax_i = add_axes_from_panel(fig, panel_i, layout_fig_w, layout_fig_h)

        cols_jk = int(np.ceil(n_genes / 2))
        gs_j = add_gridspec_from_panel(
            fig, panel_j, layout_fig_w, layout_fig_h, 2, cols_jk, wspace=0.15, hspace=0.18
        )
        gs_k = add_gridspec_from_panel(
            fig, panel_k, layout_fig_w, layout_fig_h, 2, cols_jk, wspace=0.15, hspace=0.18
        )
        axes_j = [fig.add_subplot(gs_j[i // cols_jk, i % cols_jk]) for i in range(2 * cols_jk)]
        axes_k = [fig.add_subplot(gs_k[i // cols_jk, i % cols_jk]) for i in range(2 * cols_jk)]
    else:
        if args.compact:
            height_ratios = [1.00, 1.00, 1.16, 1.06, 1.20]
            width_ratios = [2.05, 1.50, 1.10]
            outer_hspace = 0.26
            outer_wspace = 0.20
            ab_wspace = 0.08
            row5_wspace = 0.18
            jk_wspace = 0.14
            jk_hspace = 0.18
            fig = plt.figure(figsize=(args.fig_width, args.fig_height), facecolor="#EBEBEB")
            outer = GridSpec(
                nrows=5,
                ncols=3,
                figure=fig,
                height_ratios=height_ratios,
                width_ratios=width_ratios,
                hspace=outer_hspace,
                wspace=outer_wspace,
            )
            fig.subplots_adjust(left=0.045, right=0.985, top=0.985, bottom=0.045)
            slot_a = outer[0, 0:2]
            slot_c = outer[0, 2]
            slot_b = outer[1, 0:2]
            slot_d = outer[2, 0]
            slot_e = outer[2, 1]
            slot_f = outer[2, 2]
            slot_g = outer[3, 0]
            slot_h = outer[3, 1]
            slot_i = outer[3, 2]
            slot_row5 = outer[4, :]
            ab_wspace_full = ab_wspace
        else:
            # Use one unified page grid for a dense, reference-like composition.
            height_ratios = [1.35, 1.35, 1.25, 1.00, 1.05]
            width_ratios = [2.18, 1.60, 1.16]
            outer_hspace = 0.16
            outer_wspace = 0.20
            row5_wspace = 0.20
            jk_wspace = 0.15
            jk_hspace = 0.18
            fig = plt.figure(figsize=(args.fig_width, args.fig_height), facecolor="#EBEBEB")
            outer = GridSpec(
                nrows=5,
                ncols=3,
                figure=fig,
                height_ratios=height_ratios,
                width_ratios=width_ratios,
                hspace=outer_hspace,
                wspace=outer_wspace,
            )
            fig.subplots_adjust(left=0.030, right=0.992, top=0.992, bottom=0.035)
            slot_a = outer[0, 0:2]
            slot_c = outer[0, 2]
            slot_b = outer[1, 0:2]
            slot_d = outer[2, 0]
            slot_e = outer[2, 1]
            slot_f = outer[2, 2]
            slot_g = outer[3, 0]
            slot_h = outer[3, 1]
            slot_i = outer[3, 2]
            slot_row5 = outer[4, :]
            ab_wspace_full = 0.14

        gs_a = slot_a.subgridspec(1, n_genes, wspace=ab_wspace_full)
        axes_a = [fig.add_subplot(gs_a[0, i]) for i in range(n_genes)]
        ax_c = fig.add_subplot(slot_c)
        gs_b = slot_b.subgridspec(1, n_genes, wspace=ab_wspace_full)
        axes_b = [fig.add_subplot(gs_b[0, i]) for i in range(n_genes)]
        ax_d = fig.add_subplot(slot_d)
        gs_e = slot_e.subgridspec(1, 2, wspace=0.28)
        ax_e1 = fig.add_subplot(gs_e[0, 0])
        ax_e2 = fig.add_subplot(gs_e[0, 1])
        gs_f = slot_f.subgridspec(2, 2, wspace=0.36, hspace=0.42)
        ax_f11 = fig.add_subplot(gs_f[0, 0])
        ax_f12 = fig.add_subplot(gs_f[0, 1])
        ax_f21 = fig.add_subplot(gs_f[1, 0])
        ax_f22 = fig.add_subplot(gs_f[1, 1])
        gs_g = slot_g.subgridspec(1, 2, wspace=0.25)
        ax_g1 = fig.add_subplot(gs_g[0, 0])
        ax_g2 = fig.add_subplot(gs_g[0, 1])
        ax_h = fig.add_subplot(slot_h)
        ax_i = fig.add_subplot(slot_i)
        row5 = slot_row5.subgridspec(1, 2, wspace=row5_wspace, width_ratios=[1, 1])
        cols_jk = int(np.ceil(n_genes / 2))
        gs_j = row5[0, 0].subgridspec(2, cols_jk, wspace=jk_wspace, hspace=jk_hspace)
        gs_k = row5[0, 1].subgridspec(2, cols_jk, wspace=jk_wspace, hspace=jk_hspace)
        axes_j = [fig.add_subplot(gs_j[i // cols_jk, i % cols_jk]) for i in range(2 * cols_jk)]
        axes_k = [fig.add_subplot(gs_k[i // cols_jk, i % cols_jk]) for i in range(2 * cols_jk)]

    # A
    for i, g in enumerate(genes):
        weighted_density(
            axes_a[i],
            coords_a,
            expr[:, i],
            g,
            xlim=xlim_a,
            ylim=ylim_a,
            outline_xy=outline_a,
        )
    add_panel_letter(axes_a[0], "A")

    # C
    colors = sns.color_palette("tab10", n_colors=n_genes)
    for i, g in enumerate(genes):
        score = expr[:, i]
        fpr, tpr, _ = roc_curve(y_c, score)
        auc_i = auc(fpr, tpr)
        ax_c.plot(fpr, tpr, lw=2.5, color=colors[i], label=f"{g}, AUC={auc_i:.3f}")
    ax_c.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax_c.set_xlabel("1 - Specificity", fontsize=18)
    ax_c.set_ylabel("Sensitivity", fontsize=18)
    ax_c.tick_params(labelsize=13)
    ax_c.legend(loc="lower right", fontsize=12, frameon=False)
    add_panel_letter(ax_c, "C")

    # B
    for i, g in enumerate(genes):
        weighted_density(
            axes_b[i],
            coords_b,
            expr_malignant[:, i],
            g,
            xlim=xlim_b,
            ylim=ylim_b,
            outline_xy=outline_b,
        )
    add_panel_letter(axes_b[0], "B")

    # D
    sns.boxplot(
        data=fit_df,
        x="model",
        y="auc",
        order=model_order,
        hue="model",
        palette=model_palette,
        dodge=False,
        ax=ax_d,
        fliersize=2,
        linewidth=1,
        saturation=0.95,
    )
    if ax_d.legend_ is not None:
        ax_d.legend_.remove()
    sns.stripplot(
        data=fit_df,
        x="model",
        y="auc",
        order=model_order,
        ax=ax_d,
        color="#555555",
        size=2.0,
        alpha=0.55,
        jitter=0.15,
    )
    ax_d.set_title("Machine Learning model benchmark", fontsize=13, fontweight="bold")
    ax_d.set_xlabel("Learn Model", fontsize=9)
    ax_d.set_ylabel("Mean AUC", fontsize=9)
    ax_d.tick_params(axis="x", rotation=42, labelsize=8)
    ax_d.tick_params(axis="y", labelsize=8)
    ax_d.grid(axis="y", alpha=0.25, color="#CFCFCF")
    add_panel_letter(ax_d, "D")

    # E
    for name in model_order:
        fit = next(f for f in fits if f.name == name)
        c = model_palette[name]
        prec, rec, _ = precision_recall_curve(y_test, fit.test_proba)
        fpr, tpr, _ = roc_curve(y_test, fit.test_proba)
        ax_e1.plot(rec, prec, lw=1.8, color=c, label=name)
        ax_e2.plot(fpr, tpr, lw=1.8, color=c, label=name)
    ax_e1.set_xlim(0, 1)
    ax_e1.set_ylim(0.5, 1.01)
    ax_e1.set_xlabel("Recall", fontsize=9)
    ax_e1.set_ylabel("Precision", fontsize=9)
    ax_e1.tick_params(labelsize=8)
    ax_e2.set_xlim(0, 1)
    ax_e2.set_ylim(0, 1.01)
    ax_e2.set_xlabel("1 - Specificity", fontsize=9)
    ax_e2.set_ylabel("Sensitivity", fontsize=9)
    ax_e2.tick_params(labelsize=8)
    ax_e2.plot([0, 1], [0, 1], "--", color="lightgray", lw=1)
    ax_e2.legend(
        title="Learner",
        fontsize=7,
        title_fontsize=8,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )
    add_panel_letter(ax_e1, "E")

    # F
    fpr_tr, tpr_tr, _ = roc_curve(y_train, best_fit.train_proba)
    fpr_te, tpr_te, _ = roc_curve(y_test, best_fit.test_proba)
    auc_tr = auc(fpr_tr, tpr_tr)
    auc_te = auc(fpr_te, tpr_te)
    prec_tr, rec_tr, _ = precision_recall_curve(y_train, best_fit.train_proba)
    prec_te, rec_te, _ = precision_recall_curve(y_test, best_fit.test_proba)
    pr_tr = auc(rec_tr, prec_tr)
    pr_te = auc(rec_te, prec_te)
    line_c = "#7FD0D5"
    base_train = float(np.mean(y_train))
    base_test = float(np.mean(y_test))
    ax_f11.plot(fpr_tr, tpr_tr, color=line_c, lw=2)
    ax_f11.plot([0, 1], [0, 1], "--", color="lightgray", lw=1)
    ax_f11.set_title(f"development set N={len(y_train)}, AUC = {auc_tr:.2f}", fontsize=7)
    ax_f11.set_xlabel("False Positive Rate (FPR)", fontsize=7)
    ax_f11.set_ylabel("True Positive Rate (TPR)", fontsize=7)
    ax_f12.plot(rec_tr, prec_tr, color=line_c, lw=2)
    ax_f12.axhline(base_train, color="lightgray", lw=0.9, ls="--")
    ax_f12.set_title(f"development set N={len(y_train)}, PRAUC = {pr_tr:.2f}", fontsize=7)
    ax_f12.set_xlabel("Recall", fontsize=7)
    ax_f12.set_ylabel("Precision", fontsize=7)
    ax_f21.plot(fpr_te, tpr_te, color=line_c, lw=2)
    ax_f21.plot([0, 1], [0, 1], "--", color="lightgray", lw=1)
    ax_f21.set_title(f"test set N={len(y_test)}, AUC = {auc_te:.2f}", fontsize=7)
    ax_f21.set_xlabel("False Positive Rate (FPR)", fontsize=7)
    ax_f21.set_ylabel("True Positive Rate (TPR)", fontsize=7)
    ax_f22.plot(rec_te, prec_te, color=line_c, lw=2)
    ax_f22.axhline(base_test, color="lightgray", lw=0.9, ls="--")
    ax_f22.set_title(f"test set N={len(y_test)}, PRAUC = {pr_te:.2f}", fontsize=7)
    ax_f22.set_xlabel("Recall", fontsize=7)
    ax_f22.set_ylabel("Precision", fontsize=7)
    for ax in [ax_f11, ax_f12, ax_f21, ax_f22]:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.tick_params(labelsize=7)
    add_panel_letter(ax_f11, "F")

    # G
    format_conf_mat_reference(ax_g1, y_train, y_pred_train, "Confusion matrix for the train set")
    format_conf_mat_reference(ax_g2, y_test, y_pred_test, "Confusion matrix for the test set")
    add_panel_letter(ax_g1, "G")

    # H
    plot_panel_h_reference(ax_h, y_test, best_fit.test_proba)
    add_panel_letter(ax_h, "H")

    # I
    plot_panel_i_reference(ax_i, shap_vals, X_test_mm, genes, correct_test)
    add_panel_letter(ax_i, "I")

    for i, g in enumerate(genes):
        axj = axes_j[i]
        x = X_test_mm[:, i]
        yp = best_fit.test_proba
        plot_panel_j_reference(
            axj,
            x,
            yp,
            correct_test,
            g,
            show_ylabel=(i % cols_jk == 0),
            show_xlabel=(i >= cols_jk),
        )

        axk = axes_k[i]
        ys = shap_vals[:, i]
        plot_panel_k_reference(
            axk,
            x,
            ys,
            correct_test,
            g,
            show_ylabel=(i % cols_jk == 0),
            show_xlabel=(i >= cols_jk),
        )

    handles = [
        mpl.lines.Line2D([], [], color="black", marker="o", linestyle="None", markersize=4, label="FALSE"),
        mpl.lines.Line2D([], [], color="black", marker="^", linestyle="None", markersize=5, label="TRUE"),
    ]
    spare_j = list(range(n_genes, 2 * cols_jk))
    spare_k = list(range(n_genes, 2 * cols_jk))
    if spare_j:
        legend_ax_j = axes_j[spare_j[0]]
        legend_ax_j.axis("off")
        legend_ax_j.legend(handles=handles, title="correct_prediction", frameon=False, loc="center")
        spare_j = spare_j[1:]
    if spare_k:
        legend_ax_k = axes_k[spare_k[0]]
        legend_ax_k.axis("off")
        legend_ax_k.legend(handles=handles, title="correct_prediction", frameon=False, loc="center")
        spare_k = spare_k[1:]
    for i in spare_j:
        axes_j[i].axis("off")
    for i in spare_k:
        axes_k[i].axis("off")
    if n_genes >= 2 * cols_jk:
        fig.legend(handles=handles, title="correct_prediction", frameon=False, loc="lower center", ncol=2)
    add_panel_letter(axes_j[0], "J")
    add_panel_letter(axes_k[0], "K")

    # Keep non-density panels close to the reference page background.
    for ax in [ax_c, ax_d, ax_e1, ax_e2, ax_f11, ax_f12, ax_f21, ax_f22, ax_g1, ax_g2, ax_h, ax_i]:
        ax.set_facecolor("#EBEBEB")

    out_png = os.path.join(args.out_dir, "Fig6_combined_reference_like.png")
    out_pdf = os.path.join(args.out_dir, "Fig6_combined_reference_like.pdf")
    fig.savefig(out_png, dpi=args.dpi, facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[OK] Saved: {out_png}")
    print(f"[OK] Saved: {out_pdf}")
    print(f"[INFO] Panel A cells(all atlas): {adata.n_obs}")
    print(f"[INFO] Panel B cells(malignant subset): {adata_malignant.n_obs}")
    print(f"[INFO] Coordinate source A: {coord_source_a}")
    print(f"[INFO] Coordinate source B: {coord_source_b}")
    print(f"[INFO] Panel B metadata: {panel_b_meta_path}")
    if layout_meta_path:
        print(f"[INFO] Layout metadata: {layout_meta_path}")
    print(f"[INFO] Panel C label counts: malignant={int(np.sum(y_c == 1))}, non_malignant={int(np.sum(y_c == 0))}")
    cls_counts = meta_classifier["stem_group_strict"].astype(str).value_counts()
    print(
        "[INFO] D-K classifier cohort: "
        f"N={adata_classifier.n_obs}, "
        f"DTStem={int(cls_counts.get('DTStem', 0))}, "
        f"HighStem={int(cls_counts.get('HighStem', 0))}, "
        f"LowStem={int(cls_counts.get('LowStem', 0))}"
    )
    print(f"[INFO] D-K classifier metadata: {classifier_meta_path}")
    print(f"[INFO] Analysis cache: {analysis_cache_path}")
    print(f"[INFO] Analysis cache loaded: {cache_loaded}")
    print(f"[INFO] Gene file: {gene_file}")
    print(f"[INFO] Genes used ({n_genes}): {genes}")
    print(f"[INFO] Train N={len(y_train)}, Test N={len(y_test)}, Best model={best_fit.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
