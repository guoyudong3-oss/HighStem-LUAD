#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import binom, pearsonr

try:
    import scanpy as sc

    def read_h5ad(path: str):
        return sc.read_h5ad(path)
except Exception:
    import anndata as ad

    def read_h5ad(path: str):
        return ad.read_h5ad(path)

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt

from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import RFE, SelectFromModel
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Final screening (significant version) without overwriting outputs"
    )
    parser.add_argument(
        "--run-dir",
        required=True,
        help="Run output directory. Tables/logs will be written under this directory.",
    )
    parser.add_argument(
        "--adata-path",
        default="data/processed/discovery_malignant.h5ad",
        help="Input h5ad file for malignant expression matrix.",
    )
    parser.add_argument(
        "--meta-strict-path",
        default="results/reproduce_paper/r_data/metadata_strict.csv",
        help="Input strict metadata CSV (index=cell_id, columns include cytotrace_score/stem_group_strict).",
    )
    return parser.parse_args()


def assert_not_exists(path: str) -> None:
    if os.path.exists(path):
        raise FileExistsError(f"Refusing to overwrite existing file: {path}")


def next_available_path(path: str) -> str:
    if not os.path.exists(path):
        return path

    base, ext = os.path.splitext(path)
    i = 1
    while True:
        candidate = f"{base}.{i:03d}{ext}"
        if not os.path.exists(candidate):
            return candidate
        i += 1


def write_gene_list(path: str, genes: list[str]) -> None:
    pd.DataFrame({"gene": genes}).to_csv(path, index=False)


def to_float(x: object) -> float:
    return float(np.asarray(x).reshape(-1)[0])


def setup_logging(log_file: str) -> logging.Logger:
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    log_file = next_available_path(log_file)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("plasticstem")


BORUTA_HISTORY: list = []


def manual_boruta(
    logger: logging.Logger,
    X: np.ndarray,
    y: np.ndarray,
    gene_names: list[str],
    n_iterations: int = 20,
    random_state: int = 42,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Manual Boruta-style screening.

    1) Create shadow features (shuffled versions of real ones)
    2) Train RF on real + shadow
    3) A real feature scores a hit if importance > max(shadow_importance)
    4) Keep features that win significantly more than chance (binomial test)
    """

    logger.info("      Running manual Boruta logic...")

    n_features = X.shape[1]
    hits = np.zeros(n_features)
    rf = RandomForestClassifier(n_estimators=100, n_jobs=-1, random_state=random_state)

    for i in range(n_iterations):
        np.random.seed(random_state + i)
        X_shadow = np.apply_along_axis(np.random.permutation, 0, X)
        X_combined = np.hstack([X, X_shadow])

        rf.fit(X_combined, y)
        importances = rf.feature_importances_
        real_imp = importances[:n_features]
        shadow_imp = importances[n_features:]
        hits += real_imp > np.max(shadow_imp)
        BORUTA_HISTORY.append((real_imp.copy(), shadow_imp.copy()))

    p_values = np.asarray([1 - binom.cdf(hit - 1, n_iterations, 0.5) for hit in hits])
    confirmed = [gene_names[i] for i, p in enumerate(p_values) if p < 0.05]
    return confirmed, hits, p_values


def save_barplot(
    path: str,
    series: pd.Series,
    title: str,
    xlabel: str = "Score",
    top_n: int = 30,
) -> None:
    if series.empty:
        return
    s = series.sort_values(ascending=False).head(top_n)
    height = max(4.0, 0.25 * len(s) + 1.5)
    plt.figure(figsize=(8, height))
    plt.barh(list(reversed(s.index.tolist())), list(reversed(s.values.tolist())))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def save_upset(
    path: str, selected_sets: dict[str, set[str]], max_intersections: int = 15
) -> None:
    if not selected_sets:
        return

    canonical_order = [
        ("RandomForest", ("RandomForest", "RF")),
        ("DT", ("DT",)),
        ("GBM", ("GBM",)),
        ("ABESS", ("ABESS", "ABESS_Proxy", "ABESS_truth")),
        ("LASSO", ("LASSO",)),
        ("Boruta", ("Boruta",)),
    ]
    names: list[str] = []
    display_to_source: dict[str, str] = {}
    for display_name, aliases in canonical_order:
        for alias in aliases:
            if alias in selected_sets:
                names.append(display_name)
                display_to_source[display_name] = alias
                break

    if not names:
        return

    all_genes = set()
    for display_name in names:
        all_genes.update(selected_sets[display_to_source[display_name]])
    if not all_genes:
        return

    counts: Counter[tuple[bool, ...]] = Counter()
    for g in all_genes:
        mask = tuple(g in selected_sets[display_to_source[name]] for name in names)
        counts[mask] += 1

    top = sorted(
        counts.items(),
        key=lambda item: (-item[1], -sum(item[0]), item[0]),
    )[:max_intersections]
    masks = [m for m, _ in top]
    sizes = [c for _, c in top]
    x = np.arange(len(sizes))

    fig = plt.figure(figsize=(max(4.8, 0.42 * len(sizes) + 2.2), 5.2))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.85, 4.15],
        height_ratios=[2.6, 1.35],
        wspace=0.72,
        hspace=0.06,
    )
    ax_blank = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_set = fig.add_subplot(gs[1, 0])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)

    ax_blank.axis("off")

    ax_bar.bar(x, sizes, color="#333333", width=0.78)
    ax_bar.set_ylabel("Intersection Size")
    ax_bar.set_xticks([])
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.spines["bottom"].set_visible(False)
    ax_bar.grid(False)
    ymax = max(sizes) if sizes else 1
    ax_bar.set_ylim(0, ymax * 1.12)
    for xi, yi in zip(x, sizes):
        ax_bar.text(
            xi,
            yi + ymax * 0.03,
            str(int(yi)),
            ha="center",
            va="bottom",
            fontsize=8,
            color="#333333",
        )

    set_sizes = np.array(
        [len(selected_sets[display_to_source[name]]) for name in names], dtype=float
    )
    set_colors = {
        "RandomForest": "#E41A1C",
        "DT": "#F29CA3",
        "GBM": "#33A02C",
        "ABESS": "#B2DF8A",
        "LASSO": "#1F78B4",
        "Boruta": "#A6CEE3",
    }
    y_pos = np.arange(len(names))
    ax_set.barh(
        y_pos,
        set_sizes,
        color=[set_colors.get(name, "#999999") for name in names],
        height=0.78,
        edgecolor="none",
    )
    ax_set.set_xlabel("Set Size")
    ax_set.set_yticks([])
    ax_set.invert_xaxis()
    ax_set.set_xlim(float(np.max(set_sizes) * 1.08), 0)
    ax_set.set_ylim(len(names) - 0.5, -0.5)
    ax_set.spines["top"].set_visible(False)
    ax_set.spines["right"].set_visible(False)
    ax_set.spines["left"].set_visible(False)
    ax_set.grid(False)

    ax_mat.set_xlim(-0.5, len(sizes) - 0.5)
    ax_mat.set_ylim(len(names) - 0.5, -0.5)
    for i in range(len(sizes)):
        ax_mat.scatter(
            np.repeat(i, len(names)),
            y_pos,
            color="#DDDDDD",
            s=18,
            zorder=1,
        )
    for i, mask in enumerate(masks):
        ys = [j for j, v in enumerate(mask) if v]
        if ys:
            ax_mat.scatter([i] * len(ys), ys, color="#222222", s=26, zorder=3)
        if len(ys) > 1:
            ax_mat.plot([i, i], [min(ys), max(ys)], color="#222222", linewidth=1.1, zorder=2)

    ax_mat.set_yticks(y_pos)
    ax_mat.set_yticklabels(names)
    ax_mat.set_xlabel("")
    ax_mat.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_mat.tick_params(axis="y", pad=8)
    ax_mat.spines["top"].set_visible(False)
    ax_mat.spines["right"].set_visible(False)
    ax_mat.spines["bottom"].set_visible(False)
    ax_mat.grid(False)

    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.10, top=0.98, wspace=0.72, hspace=0.08)
    plt.savefig(path)
    plt.close()


def main() -> int:
    args = parse_args()
    run_dir = os.path.abspath(args.run_dir)

    tables_dir = os.path.join(run_dir, "tables")
    logs_dir = os.path.join(run_dir, "logs")
    figures_dir = os.path.join(run_dir, "figures")
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(os.environ.get("MPLCONFIGDIR", "/tmp/mpl"), exist_ok=True)

    print(f"Using run dir: {run_dir}")

    logger = setup_logging(os.path.join(logs_dir, "03_ml_screening_final.log"))

    logger.info(
        "================================================================================"
    )
    logger.info("Stage 4: Final methodology-aligned screening (5-gene consensus)")
    logger.info(
        "================================================================================"
    )

    adata_path = os.path.abspath(args.adata_path)
    if not os.path.exists(adata_path):
        logger.error("File not found: %s", adata_path)
        return 1

    logger.info("Loading data...")
    adata = read_h5ad(adata_path)

    meta_path = os.path.abspath(args.meta_strict_path)
    if not os.path.exists(meta_path):
        logger.error("Metadata strict not found: %s", meta_path)
        return 1

    logger.info("Using adata_path: %s", adata_path)
    logger.info("Using meta_strict_path: %s", meta_path)

    meta_strict = pd.read_csv(meta_path, index_col=0)
    meta_index = set(meta_strict.index.astype(str).tolist())
    common_cells = [c for c in adata.obs_names.astype(str).tolist() if c in meta_index]
    common_cells_arr = np.asarray(common_cells, dtype=str)
    adata = adata[common_cells_arr, :].copy()
    adata.obs = adata.obs.merge(
        meta_strict[["cytotrace_score", "stem_group_strict"]],
        left_index=True,
        right_index=True,
        how="left",
    )
    use_raw = True
    if use_raw:
        logger.info("Using adata.raw expression matrix for screening and ML.")
        expr_genes = set(map(str, adata.raw.var_names))
    else:
        logger.info("adata.raw unavailable; falling back to adata.X expression matrix.")
        expr_genes = set(map(str, adata.var_names))

    logger.info("Step 4.1: Loading hdWGCNA hubs...")
    hub_path = os.path.join(tables_dir, "wgcna_hub_genes.csv")
    if not os.path.exists(hub_path):
        logger.error("WGCNA hub genes not found: %s", hub_path)
        return 1

    hub_genes_df = pd.read_csv(hub_path)
    candidate_genes = [g for g in hub_genes_df["gene"].tolist() if str(g) in expr_genes]
    logger.info("   Initial candidates from hdWGCNA: %d", len(candidate_genes))

    out_candidates = os.path.join(tables_dir, "wgcna_hub_genes_candidates_in_adata.csv")
    write_gene_list(out_candidates, candidate_genes)

    logger.info("Step 4.2: Pearson correlation screening (within HighStem)...")
    high_stem_adata = adata[adata.obs["stem_group_strict"] == "HighStem"].copy()
    if use_raw:
        expr_matrix = high_stem_adata.raw[:, candidate_genes].X
    else:
        expr_matrix = high_stem_adata[:, candidate_genes].X
    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()
    expr_array = np.asarray(expr_matrix)
    scores = high_stem_adata.obs["cytotrace_score"].values

    pearson_rows: list[dict[str, object]] = []
    r_thresh = 0.15
    p_thresh = 0.05
    for i, gene in enumerate(candidate_genes):
        gene_std = float(np.std(expr_array[:, i]))
        if gene_std == 0:
            pearson_rows.append(
                {
                    "gene": gene,
                    "r": np.nan,
                    "p": np.nan,
                    "std": gene_std,
                    "pass_filter": False,
                    "reason": "std=0",
                }
            )
            continue

        r_val, p_val = pearsonr(expr_array[:, i], scores)
        r = to_float(r_val)
        p = to_float(p_val)
        pass_filter = bool(r > r_thresh and p < p_thresh)
        pearson_rows.append(
            {
                "gene": gene,
                "r": float(r),
                "p": float(p),
                "std": gene_std,
                "pass_filter": pass_filter,
                "reason": "" if pass_filter else f"r<={r_thresh} or p>={p_thresh}",
            }
        )

    pearson_df = pd.DataFrame(pearson_rows)
    out_pearson_all = os.path.join(tables_dir, "high_stem_pearson_all_candidates.csv")
    pearson_df.to_csv(out_pearson_all, index=False)
    logger.info("   Saved Pearson full table to %s", out_pearson_all)

    if "pass_filter" not in pearson_df.columns:
        logger.error("No candidate genes available for Pearson screening.")
        return 1

    activity_df = pearson_df[pearson_df["pass_filter"]].copy()
    sort_idx = np.argsort(np.asarray(activity_df["r"]))[::-1]
    activity_df = activity_df.iloc[sort_idx].copy()
    logger.info(
        "   Total positively correlated genes (p<0.05, r>0.15): %d",
        len(activity_df),
    )
    activity_genes = activity_df["gene"].tolist()
    logger.info("   HighStem activity gene set size: %d", len(activity_genes))

    out_activity = os.path.join(tables_dir, "high_stem_activity_genes_final.csv")
    activity_df[["gene", "r", "p", "std"]].to_csv(out_activity, index=False)

    out_fig5a = os.path.join(figures_dir, "Fig5A_pearson_top100.pdf")
    if not activity_df.empty:
        pearson_series = pd.Series(
            activity_df["r"].values, index=activity_df["gene"].astype(str).values
        )
        save_barplot(
            out_fig5a,
            pearson_series,
            title="Top correlated genes (Pearson r)",
            xlabel="Pearson r",
            top_n=min(100, len(pearson_series)),
        )

    logger.info("Step 4.3: Running 6 ML algorithms with intrinsic selection...")
    train_adata = adata[
        adata.obs["stem_group_strict"].isin(["HighStem", "LowStem"]), :
    ].copy()
    if use_raw:
        X = train_adata.raw[:, activity_genes].X
    else:
        X = train_adata[:, activity_genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X_array = np.asarray(X)
    y = (train_adata.obs["stem_group_strict"] == "HighStem").astype(int).values

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_array)

    selected_sets: dict[str, set[str]] = {}
    feature_scores: dict[str, pd.Series] = {}
    random_seed = 42

    logger.info("   [1/6] Training LASSO...")
    lasso = LassoCV(cv=5, max_iter=10000, random_state=random_seed).fit(X_scaled, y)
    lasso_genes = [
        activity_genes[i] for i, coef in enumerate(lasso.coef_) if abs(coef) > 1e-5
    ]
    selected_sets["LASSO"] = set(lasso_genes)
    feature_scores["LASSO"] = pd.Series(
        np.abs(lasso.coef_), index=pd.Index(activity_genes, name="gene")
    )
    logger.info("      LASSO selected %d genes.", len(lasso_genes))
    write_gene_list(
        os.path.join(tables_dir, "ml_selected_genes_LASSO.csv"), lasso_genes
    )

    logger.info("   [2/6] Training Boruta...")
    boruta_genes, boruta_hits, boruta_pvals = manual_boruta(
        logger, X_scaled, y, activity_genes, random_state=random_seed
    )
    selected_sets["Boruta"] = set(boruta_genes)
    feature_scores["Boruta"] = pd.Series(
        boruta_hits, index=pd.Index(activity_genes, name="gene")
    )
    logger.info("      Boruta selected %d genes.", len(boruta_genes))
    write_gene_list(
        os.path.join(tables_dir, "ml_selected_genes_Boruta.csv"), sorted(boruta_genes)
    )
    boruta_df = pd.DataFrame(
        {
            "gene": activity_genes,
            "hits": boruta_hits.astype(int),
            "p_value": boruta_pvals,
            "confirmed": boruta_pvals < 0.05,
        }
    )
    boruta_df.sort_values(["hits", "p_value"], ascending=[False, True]).to_csv(
        os.path.join(tables_dir, "ml_boruta_scores.csv"), index=False
    )

    logger.info("   [3/6] Training Random Forest...")
    rf = RandomForestClassifier(n_estimators=200, random_state=random_seed).fit(
        X_scaled, y
    )
    selector = SelectFromModel(rf, threshold="mean", prefit=True)
    rf_support = selector.get_support(indices=True)
    if rf_support is None:
        raise RuntimeError("RandomForest selector returned None support indices")
    rf_genes = [activity_genes[int(i)] for i in rf_support]
    selected_sets["RF"] = set(rf_genes)
    feature_scores["RF"] = pd.Series(
        rf.feature_importances_, index=pd.Index(activity_genes, name="gene")
    )
    logger.info("      RF selected %d genes.", len(rf_genes))
    write_gene_list(os.path.join(tables_dir, "ml_selected_genes_RF.csv"), rf_genes)

    logger.info("   [4/6] Training GBM...")
    gbm = GradientBoostingClassifier(random_state=random_seed).fit(X_scaled, y)
    selector = SelectFromModel(gbm, threshold="mean", prefit=True)
    gbm_support = selector.get_support(indices=True)
    if gbm_support is None:
        raise RuntimeError("GBM selector returned None support indices")
    gbm_genes = [activity_genes[int(i)] for i in gbm_support]
    selected_sets["GBM"] = set(gbm_genes)
    feature_scores["GBM"] = pd.Series(
        gbm.feature_importances_, index=pd.Index(activity_genes, name="gene")
    )
    logger.info("      GBM selected %d genes.", len(gbm_genes))
    write_gene_list(os.path.join(tables_dir, "ml_selected_genes_GBM.csv"), gbm_genes)

    logger.info("   [5/6] Training Decision Tree...")
    dt = DecisionTreeClassifier(random_state=random_seed).fit(X_scaled, y)
    selector = SelectFromModel(dt, threshold="mean", prefit=True)
    dt_support = selector.get_support(indices=True)
    if dt_support is None:
        raise RuntimeError("DecisionTree selector returned None support indices")
    dt_genes = [activity_genes[int(i)] for i in dt_support]
    selected_sets["DT"] = set(dt_genes)
    feature_scores["DT"] = pd.Series(
        dt.feature_importances_, index=pd.Index(activity_genes, name="gene")
    )
    logger.info("      DT selected %d genes.", len(dt_genes))
    write_gene_list(os.path.join(tables_dir, "ml_selected_genes_DT.csv"), dt_genes)

    logger.info("   [6/6] Training Best Subset (ABESS proxy via RFE)...")
    abess_proxy = RFE(estimator=SVC(kernel="linear"), n_features_to_select=15)
    abess_proxy.fit(X_scaled, y)
    abess_support = np.where(abess_proxy.support_)[0]
    abess_genes = [activity_genes[int(i)] for i in abess_support]
    selected_sets["ABESS_Proxy"] = set(abess_genes)
    max_rank = int(np.max(abess_proxy.ranking_)) if hasattr(abess_proxy, "ranking_") else 1
    abess_scores = (max_rank - abess_proxy.ranking_ + 1).astype(float)
    feature_scores["ABESS_Proxy"] = pd.Series(
        abess_scores, index=pd.Index(activity_genes, name="gene")
    )
    logger.info("      ABESS proxy selected %d genes.", len(abess_genes))
    write_gene_list(
        os.path.join(tables_dir, "ml_selected_genes_ABESS_Proxy.csv"), abess_genes
    )

    # ---- copy-only: export scores / LASSO path+cv / Boruta history / RFE path ----
    from sklearn.linear_model import lasso_path

    for name, ser in feature_scores.items():
        if name == "Boruta":
            continue  # ml_boruta_scores.csv (hits, p_value, confirmed) already written above; avoid case-insensitive name clash
        sc_df = ser.rename("score").reset_index()
        sc_df.columns = ["gene", "score"]
        sc_df["selected"] = sc_df["gene"].isin(selected_sets[name])
        sc_df.to_csv(os.path.join(tables_dir, f"ml_{name}_scores.csv"), index=False)
    lasso_signed = pd.DataFrame({"gene": activity_genes, "coef_signed": lasso.coef_, "coef_abs": np.abs(lasso.coef_)})
    lasso_signed.to_csv(os.path.join(tables_dir, "ml_LASSO_signed_scores.csv"), index=False)

    alphas = np.asarray(lasso.alphas_)
    _, coefs, _ = lasso_path(X_scaled, y, alphas=alphas, max_iter=10000)
    nz = (np.abs(coefs) > 1e-12).sum(axis=0)
    path_rows = []
    for li in range(len(alphas)):
        for gi, g in enumerate(activity_genes):
            path_rows.append((g, li + 1, float(alphas[li]), int(nz[li]), float(coefs[gi, li])))
    pd.DataFrame(path_rows, columns=["gene", "lambda_idx", "lambda", "nzero", "coefficient"]).to_csv(
        os.path.join(tables_dir, "ml_lasso_path.csv"), index=False
    )
    cvm = lasso.mse_path_.mean(axis=1)
    cvsd = lasso.mse_path_.std(axis=1, ddof=1) / np.sqrt(lasso.mse_path_.shape[1])
    i_min = int(np.argmin(cvm))
    ok = np.where(cvm <= cvm[i_min] + cvsd[i_min])[0]
    lambda_1se = float(alphas[ok].max())
    pd.DataFrame(
        {
            "lambda_idx": np.arange(1, len(alphas) + 1),
            "lambda": alphas,
            "cvm": cvm,
            "cvsd": cvsd,
            "nzero": nz,
            "lambda_min": float(lasso.alpha_),
            "lambda_1se": lambda_1se,
            "nobs": X_scaled.shape[0],
        }
    ).to_csv(os.path.join(tables_dir, "ml_lasso_cv.csv"), index=False)

    hist = {"iteration": np.arange(1, len(BORUTA_HISTORY) + 1)}
    for gi, g in enumerate(activity_genes):
        hist[g] = [h[0][gi] for h in BORUTA_HISTORY]
    for si in range(len(activity_genes)):
        hist[f"shadow_{si + 1}"] = [h[1][si] for h in BORUTA_HISTORY]
    pd.DataFrame(hist).to_csv(os.path.join(tables_dir, "ml_boruta_history.csv"), index=False)

    # RFE(SVC linear) elimination path, step=1 (same estimator as the proxy, recorded at every support size)
    remaining = list(range(len(activity_genes)))
    rfe_rows, rfe_tune = [], []
    while len(remaining) >= 1:
        svc = SVC(kernel="linear").fit(X_scaled[:, remaining], y)
        coef = svc.coef_.ravel()
        for k, gi in enumerate(remaining):
            rfe_rows.append((activity_genes[gi], len(remaining), float(coef[k])))
        rfe_tune.append((len(remaining), float(np.mean(svc.predict(X_scaled[:, remaining]) != y))))
        if len(remaining) == 1:
            break
        remaining.pop(int(np.argmin(np.abs(coef))))
    pd.DataFrame(rfe_rows, columns=["gene", "support_size", "coefficient"]).to_csv(
        os.path.join(tables_dir, "ml_abess_proxy_path.csv"), index=False
    )
    pd.DataFrame(rfe_tune, columns=["support_size", "tune_value"]).to_csv(
        os.path.join(tables_dir, "ml_abess_proxy_tune.csv"), index=False
    )
    # ---- end copy-only ----

    logger.info("Step 4.4: Analyzing consensus...")
    all_selected: list[str] = []
    for s in selected_sets.values():
        all_selected.extend(list(s))
    counts = Counter(all_selected)

    logger.info("   Vote counts for top genes:")
    for g, c in counts.most_common(10):
        logger.info("      %s: %d/6 models", g, c)

    vote_df = pd.DataFrame(
        [{"gene": g, "votes": int(c)} for g, c in counts.most_common()]  # sorted
    )
    out_votes = os.path.join(tables_dir, "ml_gene_vote_counts.csv")
    vote_df.to_csv(out_votes, index=False)
    logger.info("   Saved vote counts to %s", out_votes)

    # Figure outputs aligned with original Fig.5 (selection panels)
    save_barplot(
        os.path.join(figures_dir, "Fig5E_boruta_hits.pdf"),
        pd.Series(boruta_hits, index=pd.Index(activity_genes, name="gene")),
        title="Boruta hits (top features)",
        xlabel="Hit count",
        top_n=30,
    )
    boruta_p = np.clip(boruta_pvals, 1e-300, 1.0)
    save_barplot(
        os.path.join(figures_dir, "Fig5F_boruta_pvals.pdf"),
        pd.Series(-np.log10(boruta_p), index=pd.Index(activity_genes, name="gene")),
        title="Boruta significance (-log10 p)",
        xlabel="-log10(p)",
        top_n=30,
    )
    if "DT" in feature_scores:
        save_barplot(
            os.path.join(figures_dir, "Fig5G_dt_importance.pdf"),
            feature_scores["DT"],
            title="Decision Tree feature importance",
            xlabel="Importance",
            top_n=30,
        )
    if "LASSO" in feature_scores:
        save_barplot(
            os.path.join(figures_dir, "Fig5H_lasso_coefficients.pdf"),
            feature_scores["LASSO"],
            title="LASSO coefficients (abs)",
            xlabel="|Coefficient|",
            top_n=30,
        )
    if "RF" in feature_scores:
        save_barplot(
            os.path.join(figures_dir, "Fig5I_rf_importance.pdf"),
            feature_scores["RF"],
            title="Random Forest feature importance",
            xlabel="Importance",
            top_n=30,
        )
    if "GBM" in feature_scores:
        save_barplot(
            os.path.join(figures_dir, "Fig5J_gbm_importance.pdf"),
            feature_scores["GBM"],
            title="GBM feature importance",
            xlabel="Importance",
            top_n=30,
        )
    if "ABESS_Proxy" in feature_scores:
        save_barplot(
            os.path.join(figures_dir, "Fig5K_abess_proxy_rfe.pdf"),
            feature_scores["ABESS_Proxy"],
            title="ABESS proxy (RFE) ranking",
            xlabel="Score",
            top_n=30,
        )
    save_upset(
        os.path.join(figures_dir, "Fig5L_upset_intersections.pdf"), selected_sets
    )

    final_signature: list[str] = []
    for thresh in [6, 5]:
        consensus = [g for g, c in counts.items() if c >= thresh]
        if len(consensus) >= 3:
            logger.info("   Using robust consensus with >= %d algorithms.", thresh)
            final_signature = consensus
            break

    if not final_signature:
        logger.warning(
            "   Strict consensus failed. Falling back to majority vote (>=4)."
        )
        final_signature = [g for g, c in counts.items() if c >= 4]

    logger.info("   Final result: %d genes.", len(final_signature))
    logger.info("   Genes: %s", final_signature)

    expected_genes = {"S100A16", "TIMP1", "CST6", "MUC5B", "MARCKSL1"}
    if set(final_signature) == expected_genes:
        logger.info("Verification passed: matches expected 5-gene signature.")
    else:
        diff = set(final_signature) ^ expected_genes
        logger.warning("Verification differs from expected. Diff: %s", diff)

    out_final = os.path.join(tables_dir, "final_consensus_genes.csv")
    write_gene_list(out_final, final_signature)
    logger.info("Pipeline complete. Final genes saved to %s", out_final)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
