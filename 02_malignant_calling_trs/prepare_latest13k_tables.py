#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Derive newmal_stem_group_map.csv and metadata_strict_latest13k.csv from stemness_cytotrace_malignant_runfinal CSV."
    )
    p.add_argument("--stemness-csv", required=True, help="Input stemness CSV from 00d1.")
    p.add_argument("--out-map-csv", required=True, help="Output malignant map CSV.")
    p.add_argument("--out-meta-csv", required=True, help="Output strict metadata CSV.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    stemness_csv = Path(args.stemness_csv)
    out_map_csv = Path(args.out_map_csv)
    out_meta_csv = Path(args.out_meta_csv)

    if not stemness_csv.exists():
        raise FileNotFoundError(f"Missing stemness CSV: {stemness_csv}")

    df = pd.read_csv(stemness_csv)
    required = {"cell_id", "cytotrace_score", "stem_group"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Stemness CSV missing columns: {sorted(missing)}")

    df = df.loc[:, ["cell_id", "cytotrace_score", "stem_group"]].copy()
    df = df.dropna(subset=["cell_id", "cytotrace_score", "stem_group"])
    df["cell_id"] = df["cell_id"].astype(str)
    df["stem_group"] = df["stem_group"].astype(str)

    map_df = df.loc[:, ["cell_id", "stem_group"]].copy()
    meta_df = df.rename(columns={"stem_group": "stem_group_strict"}).copy()

    out_map_csv.parent.mkdir(parents=True, exist_ok=True)
    out_meta_csv.parent.mkdir(parents=True, exist_ok=True)
    map_df.to_csv(out_map_csv, index=False)
    meta_df.to_csv(out_meta_csv, index=False)

    print(f"OUT_MAP={out_map_csv.resolve()}")
    print(f"OUT_META={out_meta_csv.resolve()}")
    print(f"N_ROWS={len(df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
