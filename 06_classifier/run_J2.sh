#!/bin/bash
cd "${PROJECT_ROOT:?set PROJECT_ROOT to the repository root}"
J2=results/reproduce_paper/runs/final_malignant_v3_20261009/J2; P=06_classifier; PB=results/reproduce_paper/runs/paper_main_epiref_latest13k_noclobber_flow_20260213_152245/tables/cnv_paper_exact_cell_metadata_runfinal_ref_mean_3sd.csv; L=06_classifier/run_with_sklearn_labels.py; PY="${PY_BIN:-python3}"
($PY $L $P/make_fig6_reference_like.py --run-dir $J2/run --out-dir $J2/figures_full --meta-path-panel-b $PB --meta-path-classifier $J2/run/tables/metadata_strict_latest13k.csv > $J2/full.log 2>&1; echo "full $?" >> $J2/status.txt) &
($PY $L $P/make_fig6_panel_d_only.py --run-dir $J2/run --meta-path-classifier $J2/run/tables/metadata_strict_latest13k.csv --out-dir $J2/panel_D > $J2/panel_D.log 2>&1; echo "D $?" >> $J2/status.txt) &
($PY $L $P/make_fig6_panel_d_only.py --run-dir $J2/run_regression --meta-path-classifier $J2/run_regression/tables/metadata_strict_latest13k.csv --out-dir $J2/panel_D_regression_origorder > $J2/panel_D_regression.log 2>&1; echo "Dreg $?" >> $J2/status.txt) &
wait
