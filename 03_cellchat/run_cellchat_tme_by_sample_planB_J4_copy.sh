#!/usr/bin/env bash
set -euo pipefail

# Plan B:
# 1) Per-sample rerun of CellChat (same labeling logic as full-TME run)
# 2) Sample-level statistical tests on hub/HighStem axes

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${PROJECT_ROOT:?set PROJECT_ROOT to the repository root}"
cd "${ROOT_DIR}"

RUN_NAME="${1:-paper_main_tme_cellchat_by_sample_planB_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="results/reproduce_paper/runs/${RUN_NAME}"

SOURCE_RUN="${SOURCE_RUN:-results/reproduce_paper/runs/final_malignant_v3_20261009/J4/main_1200}"
SOURCE_META_FULL="${SOURCE_META_FULL:-${SOURCE_RUN}/cellchat_input/cellchat_metadata_full.csv}"
SOURCE_MAP_CSV="${SOURCE_MAP_CSV:-results/reproduce_paper/runs/paper_main_tme_cellchat_highratio10_20260214_095043/tables/newmal_stem_group_map.csv}"

INTEGRATED_H5AD="${INTEGRATED_H5AD:-data/processed/discovery_integrated.h5ad}"
MALIGNANT_RDS="${MALIGNANT_RDS:-results/reproduce_paper/runs/paper_main_epiref_full_noplot_20260209_114337/rds/seurat_malignant_with_stemness_runfinal_ref_mean_3sd.rds}"

PY_BIN="${PY_BIN:-.venv311/bin/python}"
SAMPLE_COL="${SAMPLE_COL:-sample_id}"
LABEL_MODE="${LABEL_MODE:-original}"
CELLTYPE_ANNOT_COL="${CELLTYPE_ANNOT_COL:-celltype_annot}"

# Sample selection thresholds (eligibility to rerun sample-level CellChat)
MIN_HIGH_CELLS="${MIN_HIGH_CELLS:-5}"
MIN_ENDO_CELLS="${MIN_ENDO_CELLS:-5}"
MIN_TOTAL_CELLS="${MIN_TOTAL_CELLS:-200}"
MIN_GROUPS_PRESENT="${MIN_GROUPS_PRESENT:-4}"
SAMPLE_LIMIT="${SAMPLE_LIMIT:-0}"   # 0 means no limit

# 18a thresholds within each sample
MAX_PER_GROUP="${MAX_PER_GROUP:-0}"              # keep all cells within sample
MIN_CELLS_PER_GROUP="${MIN_CELLS_PER_GROUP:-5}"  # drop tiny groups within sample

if [[ -e "${RUN_DIR}" ]]; then
  echo "[ERROR] Run directory already exists: ${RUN_DIR}" >&2
  exit 1
fi

for f in "${SOURCE_META_FULL}" "${SOURCE_MAP_CSV}" "${INTEGRATED_H5AD}" "${MALIGNANT_RDS}"; do
  if [[ ! -f "${f}" ]]; then
    echo "[ERROR] Missing required input: ${f}" >&2
    exit 1
  fi
done

if [[ ! -x "${PY_BIN}" ]]; then
  echo "[ERROR] Python executable not found or not executable: ${PY_BIN}" >&2
  exit 1
fi

mkdir -p "${RUN_DIR}/tables" "${RUN_DIR}/logs" "${RUN_DIR}/by_sample"
cp "${SOURCE_MAP_CSV}" "${RUN_DIR}/tables/newmal_stem_group_map.csv"

SAMPLE_COUNTS_CSV="${RUN_DIR}/tables/sample_selection_counts.csv"
ELIGIBLE_LIST_TXT="${RUN_DIR}/tables/sample_selection_eligible.txt"
SUCCESS_LIST_TXT="${RUN_DIR}/tables/sample_runs_success.txt"
FAILED_LIST_TXT="${RUN_DIR}/tables/sample_runs_failed.txt"
: > "${SUCCESS_LIST_TXT}"
: > "${FAILED_LIST_TXT}"

echo "[INFO] Run dir: ${RUN_DIR}"
echo "[INFO] Source metadata: ${SOURCE_META_FULL}"
echo "[INFO] Selecting eligible samples..."
Rscript -e '
args <- commandArgs(trailingOnly = TRUE)
meta_file <- args[1]
out_counts <- args[2]
out_list <- args[3]
min_high <- as.integer(args[4])
min_endo <- as.integer(args[5])
min_total <- as.integer(args[6])
min_groups <- as.integer(args[7])
sample_limit <- as.integer(args[8])

library(data.table)
dt <- fread(meta_file)
need <- c("sample_id", "cellchat_label")
if (!all(need %in% colnames(dt))) {
  stop("metadata_full must include columns: sample_id, cellchat_label")
}
dt <- dt[!is.na(sample_id) & !is.na(cellchat_label)]
cnt <- dt[, .N, by = .(sample_id, cellchat_label)]
wide <- dcast(cnt, sample_id ~ cellchat_label, value.var = "N", fill = 0)

if (!("HighStem" %in% colnames(wide))) wide[, HighStem := 0L]
if (!("Endothelial cells" %in% colnames(wide))) wide[, `Endothelial cells` := 0L]

label_cols <- setdiff(colnames(wide), "sample_id")
wide[, total_cells := rowSums(.SD), .SDcols = label_cols]
wide[, n_groups_present := rowSums(.SD > 0), .SDcols = label_cols]

elig <- wide[
  HighStem >= min_high &
  `Endothelial cells` >= min_endo &
  total_cells >= min_total &
  n_groups_present >= min_groups
]
setorder(elig, -HighStem, -`Endothelial cells`)
if (sample_limit > 0 && nrow(elig) > sample_limit) {
  elig <- elig[1:sample_limit]
}

fwrite(wide, out_counts)
writeLines(elig$sample_id, out_list)
cat("Eligible samples:", nrow(elig), "\n")
if (nrow(elig) > 0) {
  print(elig[, .(sample_id, HighStem, `Endothelial cells`, total_cells, n_groups_present)])
}
' "${SOURCE_META_FULL}" "${SAMPLE_COUNTS_CSV}" "${ELIGIBLE_LIST_TXT}" \
  "${MIN_HIGH_CELLS}" "${MIN_ENDO_CELLS}" "${MIN_TOTAL_CELLS}" "${MIN_GROUPS_PRESENT}" "${SAMPLE_LIMIT}" \
  | tee "${RUN_DIR}/logs/00_sample_selection.log"

if [[ ! -s "${ELIGIBLE_LIST_TXT}" ]]; then
  echo "[ERROR] No eligible samples after filtering. Adjust thresholds and rerun." >&2
  exit 1
fi

echo "[INFO] Starting per-sample CellChat reruns..."
while IFS= read -r SAMPLE_ID; do
  [[ -z "${SAMPLE_ID}" ]] && continue
  SAMPLE_DIR="${RUN_DIR}/by_sample/${SAMPLE_ID}"
  mkdir -p "${SAMPLE_DIR}/logs" "${SAMPLE_DIR}/tables"
  echo "------------------------------------------------------------"
  echo "[SAMPLE] ${SAMPLE_ID}"

  if "${PY_BIN}" 03_cellchat/18a_prepare_fig3_cellchat_inputs.py \
      --run-dir "${SAMPLE_DIR}" \
      --integrated-h5ad "${INTEGRATED_H5AD}" \
      --malignant-map-csv "${RUN_DIR}/tables/newmal_stem_group_map.csv" \
      --max-per-group "${MAX_PER_GROUP}" \
      --min-cells-per-group "${MIN_CELLS_PER_GROUP}" \
      --celltype-annot-col "${CELLTYPE_ANNOT_COL}" \
      --nonmalignant-label-mode "${LABEL_MODE}" \
      --expr-source raw_log1p --raw-gene-mode all_raw --raw-scale-factor 10000 --seed 42 \
      --sample-id "${SAMPLE_ID}" \
      --sample-col "${SAMPLE_COL}" \
      > "${SAMPLE_DIR}/logs/18a_prep.log" 2>&1; then
    :
  else
    echo "[WARN] 18a failed for sample ${SAMPLE_ID}"
    echo "${SAMPLE_ID}" >> "${FAILED_LIST_TXT}"
    continue
  fi

  if Rscript 03_cellchat/18b_fig3_cellchat_pathway_metabolism.R \
      --run-dir "${SAMPLE_DIR}" \
      --malignant-rds "${MALIGNANT_RDS}" \
      --skip-hi \
      > "${SAMPLE_DIR}/logs/18b_cellchat.log" 2>&1; then
    echo "${SAMPLE_ID}" >> "${SUCCESS_LIST_TXT}"
  else
    echo "[WARN] 18b failed for sample ${SAMPLE_ID}"
    echo "${SAMPLE_ID}" >> "${FAILED_LIST_TXT}"
    continue
  fi
done < "${ELIGIBLE_LIST_TXT}"

if [[ ! -s "${SUCCESS_LIST_TXT}" ]]; then
  echo "[ERROR] No successful sample reruns. Check per-sample logs." >&2
  exit 1
fi

echo "[INFO] Running sample-level statistics (18c)..."
"${PY_BIN}" 03_cellchat/18c_cellchat_sample_level_stats.py \
  --parent-run-dir "${RUN_DIR}" \
  --samples-dir "${RUN_DIR}/by_sample" \
  --output-prefix "sampleB" \
  | tee "${RUN_DIR}/logs/18c_sample_stats.log"

echo "[DONE] Plan B completed."
echo "Run dir: ${RUN_DIR}"
echo "Key outputs:"
echo "  - ${RUN_DIR}/tables/sampleB_per_sample_metrics.csv"
echo "  - ${RUN_DIR}/tables/sampleB_summary_metrics.csv"
echo "  - ${RUN_DIR}/tables/sampleB_summary_tests.csv"
echo "  - ${RUN_DIR}/tables/sampleB_summary_report.txt"
echo "  - ${RUN_DIR}/tables/sample_runs_success.txt"
echo "  - ${RUN_DIR}/tables/sample_runs_failed.txt"
