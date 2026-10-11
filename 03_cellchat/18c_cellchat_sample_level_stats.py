#!/usr/bin/env python3
"""Sample-level CellChat summary and statistical tests for Plan B."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Aggregate per-sample CellChat outputs and run sample-level tests."
    )
    p.add_argument(
        "--parent-run-dir",
        required=True,
        help="Parent run directory containing by_sample/<sample>/tables/*.csv",
    )
    p.add_argument(
        "--samples-dir",
        default="",
        help="Optional explicit samples directory (default: <parent-run-dir>/by_sample).",
    )
    p.add_argument(
        "--output-prefix",
        default="sampleB",
        help="Prefix for output files in <parent-run-dir>/tables.",
    )
    p.add_argument("--highstem-label", default="HighStem")
    p.add_argument("--lowstem-label", default="LowStem")
    p.add_argument("--dtstem-label", default="DTStem")
    p.add_argument("--endo-label", default="Endothelial cells")
    return p.parse_args()


def read_square_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.shape[0] == 0 or df.shape[1] == 0:
        return df
    # Align index/columns if needed.
    common = [x for x in df.index.astype(str) if x in set(df.columns.astype(str))]
    if len(common) > 0:
        df = df.loc[common, common]
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df.astype(float)


def rank_of_label(values: pd.Series, label: str) -> float:
    if label not in values.index:
        return np.nan
    ordered = values.sort_values(ascending=False)
    for i, g in enumerate(ordered.index.tolist(), start=1):
        if g == label:
            return float(i)
    return np.nan


def safe_wilcoxon_1samp(x: np.ndarray, alternative: str) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = x[np.abs(x) > 1e-15]
    n = int(x.size)
    if n < 3:
        return (np.nan, np.nan, n)
    try:
        res = stats.wilcoxon(x, alternative=alternative, zero_method="wilcox")
        return (float(res.statistic), float(res.pvalue), n)
    except Exception:
        return (np.nan, np.nan, n)


def safe_wilcoxon_paired(
    x: Iterable[float], y: Iterable[float], alternative: str
) -> tuple[float, float, int]:
    x = np.asarray(list(x), dtype=float)
    y = np.asarray(list(y), dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    d = x - y
    d = d[np.abs(d) > 1e-15]
    n = int(d.size)
    if n < 3:
        return (np.nan, np.nan, n)
    try:
        res = stats.wilcoxon(d, alternative=alternative, zero_method="wilcox")
        return (float(res.statistic), float(res.pvalue), n)
    except Exception:
        return (np.nan, np.nan, n)


def safe_binom_greater(success: int, total: int) -> float:
    if total <= 0:
        return np.nan
    return float(stats.binomtest(k=success, n=total, p=0.5, alternative="greater").pvalue)


def normalize_text(x: object) -> str:
    if x is None:
        return ""
    s = str(x).strip()
    if s.lower() in {"nan", "none", ""}:
        return ""
    return s


def read_sample_meta(sample_dir: Path, sample_id: str) -> dict:
    out = {"sample_id": sample_id, "patient_id": "", "sample_origin": ""}
    meta_file = sample_dir / "cellchat_input" / "cellchat_metadata_full.csv"
    if not meta_file.exists():
        return out
    try:
        head = pd.read_csv(meta_file, nrows=0)
    except Exception:
        return out
    cols = list(head.columns)
    keep = [c for c in ["sample_id", "patient_id", "Sample_Origin", "Sample"] if c in cols]
    if len(keep) == 0:
        return out
    try:
        meta = pd.read_csv(meta_file, usecols=keep)
    except Exception:
        return out

    if "sample_id" in meta.columns:
        sid = meta["sample_id"].astype(str)
        sub = meta.loc[sid == str(sample_id)].copy()
        if sub.shape[0] > 0:
            meta = sub

    if "patient_id" in meta.columns:
        vals = [normalize_text(x) for x in meta["patient_id"].tolist()]
        vals = [x for x in vals if x]
        if vals:
            out["patient_id"] = vals[0]

    for col in ["Sample_Origin", "Sample"]:
        if col in meta.columns:
            vals = [normalize_text(x) for x in meta[col].tolist()]
            vals = [x for x in vals if x]
            if vals:
                out["sample_origin"] = vals[0]
                break
    return out


def build_summary_metrics(df: pd.DataFrame) -> pd.DataFrame:
    n_runs = int(df.shape[0])
    m_hs_endo = df["hs_to_endo_weight"].notna()
    n_hs_endo = int(m_hs_endo.sum())
    m_endo_rank = df["endo_rank_in_high_outgoing"].notna()
    n_rank = int(m_endo_rank.sum())
    endo_top1_n = int(df.loc[m_endo_rank, "endo_top1_in_high_outgoing"].astype(int).sum())
    endo_top2_n = int(df.loc[m_endo_rank, "endo_top2_in_high_outgoing"].astype(int).sum())
    return pd.DataFrame(
        [
            {"metric": "n_sample_runs", "value": n_runs},
            {"metric": "n_samples_with_hs_to_endo", "value": n_hs_endo},
            {"metric": "n_samples_with_endo_rank_in_high_outgoing", "value": n_rank},
            {
                "metric": "freq_endo_top1_in_high_outgoing",
                "value": (endo_top1_n / n_rank) if n_rank > 0 else np.nan,
            },
            {
                "metric": "freq_endo_top2_in_high_outgoing",
                "value": (endo_top2_n / n_rank) if n_rank > 0 else np.nan,
            },
            {
                "metric": "freq_endo_global_hub_outgoing",
                "value": df["is_endo_global_hub_outgoing"].astype(int).mean(),
            },
            {
                "metric": "freq_high_global_hub_outgoing",
                "value": df["is_high_global_hub_outgoing"].astype(int).mean(),
            },
        ]
    )


def build_summary_tests(df: pd.DataFrame) -> pd.DataFrame:
    tests: list[dict] = []

    # Test A: Is hs->endo larger than high->other median within sample?
    diff = df["hs_to_endo_minus_other_median"].to_numpy(dtype=float)
    stat, pval, n = safe_wilcoxon_1samp(diff, alternative="greater")
    diff_nz = np.asarray(diff[np.isfinite(diff)])
    diff_nz = diff_nz[np.abs(diff_nz) > 1e-15]
    pos_n = int(np.sum(diff_nz > 0))
    sign_p = safe_binom_greater(pos_n, int(diff_nz.size))
    tests.append(
        {
            "test": "hs_to_endo_minus_other_median > 0 (wilcoxon_1samp)",
            "n": n,
            "statistic": stat,
            "pvalue": pval,
            "alternative": "greater",
            "extra": "",
        }
    )
    tests.append(
        {
            "test": "hs_to_endo_minus_other_median > 0 (sign/binom)",
            "n": int(diff_nz.size),
            "statistic": pos_n,
            "pvalue": sign_p,
            "alternative": "greater",
            "extra": "statistic=number_of_positive_samples",
        }
    )

    # Test B: HighStem vs LowStem outgoing strength across samples.
    stat, pval, n = safe_wilcoxon_paired(
        df["high_out_strength"], df["low_out_strength"], alternative="greater"
    )
    d = (
        df["high_out_strength"].to_numpy(dtype=float)
        - df["low_out_strength"].to_numpy(dtype=float)
    )
    d = d[np.isfinite(d)]
    d = d[np.abs(d) > 1e-15]
    pos_n = int(np.sum(d > 0))
    sign_p = safe_binom_greater(pos_n, int(d.size))
    tests.append(
        {
            "test": "high_out_strength > low_out_strength (wilcoxon_paired)",
            "n": n,
            "statistic": stat,
            "pvalue": pval,
            "alternative": "greater",
            "extra": "",
        }
    )
    tests.append(
        {
            "test": "high_out_strength > low_out_strength (sign/binom)",
            "n": int(d.size),
            "statistic": pos_n,
            "pvalue": sign_p,
            "alternative": "greater",
            "extra": "statistic=number_of_positive_samples",
        }
    )

    # Test C: HighStem vs DTStem outgoing strength across samples.
    stat, pval, n = safe_wilcoxon_paired(
        df["high_out_strength"], df["dt_out_strength"], alternative="greater"
    )
    d = (
        df["high_out_strength"].to_numpy(dtype=float)
        - df["dt_out_strength"].to_numpy(dtype=float)
    )
    d = d[np.isfinite(d)]
    d = d[np.abs(d) > 1e-15]
    pos_n = int(np.sum(d > 0))
    sign_p = safe_binom_greater(pos_n, int(d.size))
    tests.append(
        {
            "test": "high_out_strength > dt_out_strength (wilcoxon_paired)",
            "n": n,
            "statistic": stat,
            "pvalue": pval,
            "alternative": "greater",
            "extra": "",
        }
    )
    tests.append(
        {
            "test": "high_out_strength > dt_out_strength (sign/binom)",
            "n": int(d.size),
            "statistic": pos_n,
            "pvalue": sign_p,
            "alternative": "greater",
            "extra": "statistic=number_of_positive_samples",
        }
    )
    return pd.DataFrame(tests)


def main() -> int:
    args = parse_args()
    parent = Path(args.parent_run_dir).resolve()
    samples_dir = (
        Path(args.samples_dir).resolve() if args.samples_dir else parent / "by_sample"
    )
    tables_dir = parent / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    if not samples_dir.exists():
        raise FileNotFoundError(f"Samples directory not found: {samples_dir}")

    rows: list[dict] = []
    sample_dirs = sorted([p for p in samples_dir.iterdir() if p.is_dir()])
    for sdir in sample_dirs:
        sample_id = sdir.name
        tab = sdir / "tables"
        w_file = tab / "fig3_cellchat_weight_matrix.csv"
        c_file = tab / "fig3_cellchat_count_matrix.csv"
        if not w_file.exists() or not c_file.exists():
            continue

        w = read_square_matrix(w_file)
        c = read_square_matrix(c_file)
        if w.empty:
            continue

        groups = w.index.tolist()
        out_strength = w.sum(axis=1)
        in_strength = w.sum(axis=0)
        total_strength = out_strength + in_strength

        hub_out = out_strength.idxmax()
        hub_in = in_strength.idxmax()
        hub_total = total_strength.idxmax()

        high = args.highstem_label
        low = args.lowstem_label
        dt = args.dtstem_label
        endo = args.endo_label

        high_present = high in groups
        low_present = low in groups
        dt_present = dt in groups
        endo_present = endo in groups

        hs_out = float(out_strength.get(high, np.nan))
        hs_in = float(in_strength.get(high, np.nan))
        low_out = float(out_strength.get(low, np.nan))
        dt_out = float(out_strength.get(dt, np.nan))

        hs_to_endo_w = (
            float(w.loc[high, endo]) if (high_present and endo_present) else np.nan
        )
        hs_to_endo_c = (
            float(c.loc[high, endo]) if (high_present and endo_present and not c.empty) else np.nan
        )

        hs_out_targets = pd.Series(dtype=float)
        if high_present:
            hs_out_targets = w.loc[high].copy()
            if high in hs_out_targets.index:
                hs_out_targets = hs_out_targets.drop(index=high)

        hs_out_target_n = int(hs_out_targets.shape[0]) if high_present else 0
        hs_to_other_median = float(np.nan)
        hs_to_other_mean = float(np.nan)
        endo_rank_in_hs_out = float(np.nan)
        endo_top1_hs_out = False
        endo_top2_hs_out = False

        if hs_out_target_n > 0:
            hs_to_other_median = float(np.nanmedian(hs_out_targets.values))
            hs_to_other_mean = float(np.nanmean(hs_out_targets.values))
            endo_rank_in_hs_out = rank_of_label(hs_out_targets, endo)
            endo_top1_hs_out = bool(math.isfinite(endo_rank_in_hs_out) and endo_rank_in_hs_out <= 1.0)
            endo_top2_hs_out = bool(math.isfinite(endo_rank_in_hs_out) and endo_rank_in_hs_out <= 2.0)

        mal_groups = [g for g in [high, low, dt] if g in groups]
        high_mal_out = np.nan
        low_mal_out = np.nan
        dt_mal_out = np.nan
        if len(mal_groups) >= 2:
            w_mal = w.loc[mal_groups, mal_groups]
            mal_out = w_mal.sum(axis=1)
            high_mal_out = float(mal_out.get(high, np.nan))
            low_mal_out = float(mal_out.get(low, np.nan))
            dt_mal_out = float(mal_out.get(dt, np.nan))

        sample_meta = read_sample_meta(sdir, sample_id)
        row = {
            "sample_id": sample_id,
            "patient_id": sample_meta.get("patient_id", ""),
            "sample_origin": sample_meta.get("sample_origin", ""),
            "n_groups": int(len(groups)),
            "groups": ",".join(groups),
            "global_hub_outgoing": hub_out,
            "global_hub_incoming": hub_in,
            "global_hub_total": hub_total,
            "is_endo_global_hub_outgoing": bool(hub_out == endo),
            "is_endo_global_hub_incoming": bool(hub_in == endo),
            "is_high_global_hub_outgoing": bool(hub_out == high),
            "is_high_global_hub_incoming": bool(hub_in == high),
            "high_out_strength": hs_out,
            "high_in_strength": hs_in,
            "low_out_strength": low_out,
            "dt_out_strength": dt_out,
            "high_rank_outgoing": rank_of_label(out_strength, high),
            "high_rank_incoming": rank_of_label(in_strength, high),
            "high_rank_total": rank_of_label(total_strength, high),
            "hs_to_endo_weight": hs_to_endo_w,
            "hs_to_endo_count": hs_to_endo_c,
            "hs_to_other_median_weight": hs_to_other_median,
            "hs_to_other_mean_weight": hs_to_other_mean,
            "hs_to_endo_minus_other_median": (
                hs_to_endo_w - hs_to_other_median
                if (math.isfinite(hs_to_endo_w) and math.isfinite(hs_to_other_median))
                else np.nan
            ),
            "hs_out_target_n": hs_out_target_n,
            "endo_rank_in_high_outgoing": endo_rank_in_hs_out,
            "endo_top1_in_high_outgoing": endo_top1_hs_out,
            "endo_top2_in_high_outgoing": endo_top2_hs_out,
            "high_malignant_out_strength": high_mal_out,
            "low_malignant_out_strength": low_mal_out,
            "dt_malignant_out_strength": dt_mal_out,
        }
        rows.append(row)

    per_sample = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    if per_sample.empty:
        raise RuntimeError("No valid per-sample CellChat outputs found.")

    out_per_sample = tables_dir / f"{args.output_prefix}_per_sample_metrics.csv"
    per_sample.to_csv(out_per_sample, index=False)

    # Global summary metrics and tests.
    summary_metrics = build_summary_metrics(per_sample)
    out_summary_metrics = tables_dir / f"{args.output_prefix}_summary_metrics.csv"
    summary_metrics.to_csv(out_summary_metrics, index=False)
    tests_df = build_summary_tests(per_sample)
    out_tests = tables_dir / f"{args.output_prefix}_summary_tests.csv"
    tests_df.to_csv(out_tests, index=False)

    # Stratified analyses.
    strat_metrics_all: list[pd.DataFrame] = []
    strat_tests_all: list[pd.DataFrame] = []
    for strat_col in ["patient_id", "sample_origin"]:
        if strat_col not in per_sample.columns:
            continue
        sub_all = per_sample.loc[per_sample[strat_col].astype(str).str.len() > 0].copy()
        if sub_all.empty:
            continue

        mrows: list[pd.DataFrame] = []
        trows: list[pd.DataFrame] = []
        for sval, sdf in sub_all.groupby(strat_col):
            sm = build_summary_metrics(sdf).copy()
            sm.insert(0, "stratum", str(sval))
            sm.insert(0, "strat_col", strat_col)
            mrows.append(sm)

            st = build_summary_tests(sdf).copy()
            st.insert(0, "stratum", str(sval))
            st.insert(0, "strat_col", strat_col)
            st.insert(0, "n_sample_runs", int(sdf.shape[0]))
            trows.append(st)

        if mrows:
            out_m = pd.concat(mrows, ignore_index=True)
            out_m_file = tables_dir / f"{args.output_prefix}_stratified_metrics_by_{strat_col}.csv"
            out_m.to_csv(out_m_file, index=False)
            strat_metrics_all.append(out_m)
        if trows:
            out_t = pd.concat(trows, ignore_index=True)
            out_t_file = tables_dir / f"{args.output_prefix}_stratified_tests_by_{strat_col}.csv"
            out_t.to_csv(out_t_file, index=False)
            strat_tests_all.append(out_t)

    # Human-readable report.
    n_runs = int(per_sample.shape[0])
    m_endo_rank = per_sample["endo_rank_in_high_outgoing"].notna()
    n_rank = int(m_endo_rank.sum())
    endo_top1_n = int(per_sample.loc[m_endo_rank, "endo_top1_in_high_outgoing"].astype(int).sum())
    endo_top2_n = int(per_sample.loc[m_endo_rank, "endo_top2_in_high_outgoing"].astype(int).sum())

    top_hs_endo = (
        per_sample.loc[per_sample["hs_to_endo_weight"].notna(), ["sample_id", "hs_to_endo_weight"]]
        .sort_values("hs_to_endo_weight", ascending=False)
        .head(10)
    )
    report = tables_dir / f"{args.output_prefix}_summary_report.txt"
    with report.open("w", encoding="utf-8") as f:
        f.write("Plan B sample-level CellChat summary\n")
        f.write(f"parent_run_dir: {parent}\n")
        f.write(f"samples_dir: {samples_dir}\n")
        f.write(f"sample_runs_count: {n_runs}\n")
        f.write(f"endo_top1_in_high_outgoing: {endo_top1_n}/{n_rank}\n")
        f.write(f"endo_top2_in_high_outgoing: {endo_top2_n}/{n_rank}\n")
        f.write(
            "endo_global_hub_outgoing_freq: "
            f"{per_sample['is_endo_global_hub_outgoing'].astype(int).mean():.4f}\n"
        )
        f.write(
            "high_global_hub_outgoing_freq: "
            f"{per_sample['is_high_global_hub_outgoing'].astype(int).mean():.4f}\n"
        )
        f.write("\nTop samples by hs_to_endo_weight:\n")
        if top_hs_endo.empty:
            f.write("NA\n")
        else:
            f.write(top_hs_endo.to_string(index=False))
            f.write("\n")
        f.write("\nSummary tests:\n")
        f.write(tests_df.to_string(index=False))
        f.write("\n")
        if strat_metrics_all:
            f.write("\nStratified metrics files:\n")
            for strat_col in ["patient_id", "sample_origin"]:
                p = tables_dir / f"{args.output_prefix}_stratified_metrics_by_{strat_col}.csv"
                if p.exists():
                    f.write(f"- {p}\n")
        if strat_tests_all:
            f.write("\nStratified tests files:\n")
            for strat_col in ["patient_id", "sample_origin"]:
                p = tables_dir / f"{args.output_prefix}_stratified_tests_by_{strat_col}.csv"
                if p.exists():
                    f.write(f"- {p}\n")

    print("Saved per-sample metrics:", out_per_sample)
    print("Saved summary metrics:", out_summary_metrics)
    print("Saved summary tests:", out_tests)
    for strat_col in ["patient_id", "sample_origin"]:
        p1 = tables_dir / f"{args.output_prefix}_stratified_metrics_by_{strat_col}.csv"
        p2 = tables_dir / f"{args.output_prefix}_stratified_tests_by_{strat_col}.csv"
        if p1.exists():
            print(f"Saved stratified metrics ({strat_col}):", p1)
        if p2.exists():
            print(f"Saved stratified tests ({strat_col}):", p2)
    print("Saved text report:", report)
    print("\nTop-line:")
    print(f" - sample_runs_count = {n_runs}")
    print(f" - endo_top1_in_high_outgoing = {endo_top1_n}/{n_rank}")
    print(
        " - freq_endo_global_hub_outgoing = "
        f"{per_sample['is_endo_global_hub_outgoing'].astype(int).mean():.3f}"
    )
    print(
        " - freq_high_global_hub_outgoing = "
        f"{per_sample['is_high_global_hub_outgoing'].astype(int).mean():.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
