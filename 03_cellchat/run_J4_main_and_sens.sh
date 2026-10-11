#!/bin/bash
# CellChat main analysis (max 1,200 cells per group) and two sensitivity analyses; outputs under results/.../J4/
cd "${PROJECT_ROOT:?set PROJECT_ROOT to the repository root}"
V=results/reproduce_paper/runs/final_malignant_v3_20261009/J4; S="03_cellchat"; MAP=results/reproduce_paper/runs/paper_main_tme_cellchat_highratio10_20260214_095043/tables/newmal_stem_group_map.csv; PY="${PY_BIN:-python3}"
H5=data/processed/discovery_integrated.h5ad
COMMON="--integrated-h5ad $H5 --malignant-map-csv $MAP --min-cells-per-group 50 --nonmalignant-label-mode original --expr-source raw_log1p --raw-gene-mode all_raw --raw-scale-factor 10000 --seed 42"
t0=$(date +%s)
step(){ echo "[$(date +%H:%M:%S)] $1 (elapsed $(( ($(date +%s)-t0)/60 )) min)" >> $V/progress.log; }
step "main 18a start"
$PY "$S/18a_prepare_fig3_cellchat_inputs.py" --run-dir $V/main_1200 $COMMON --max-per-group 1200 > $V/main_1200_18a.log 2>&1 || { step "main 18a FAILED"; exit 1; }
step "main 18b start"
Rscript "$S/18b_fig3_cellchat_pathway_metabolism.R" --run-dir $V/main_1200 > $V/main_1200_18b.log 2>&1 || step "main 18b exit nonzero"
step "sensB 18a start"
$PY $S/18a_copy_fullstem.py --run-dir $V/sensB_fullstem $COMMON --max-per-group 1200 > $V/sensB_18a.log 2>&1 || { step "sensB 18a FAILED"; exit 1; }
step "sensB 18b start"
Rscript "$S/18b_fig3_cellchat_pathway_metabolism.R" --run-dir $V/sensB_fullstem --skip-hi > $V/sensB_18b.log 2>&1 || step "sensB 18b exit nonzero"
step "sensC prep (copy main inputs)"
mkdir -p $V/sensC_popsizeFALSE; cp -R $V/main_1200/cellchat_input $V/sensC_popsizeFALSE/cellchat_input
step "sensC 18b start"
Rscript $S/18b_copy_popsizeFALSE.R --run-dir $V/sensC_popsizeFALSE --skip-hi > $V/sensC_18b.log 2>&1 || step "sensC 18b exit nonzero"
step "ALL DONE"
