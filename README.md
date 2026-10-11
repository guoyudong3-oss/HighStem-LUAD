# HighStem-LUAD

Analysis code for the study *"Single-cell and spatial transcriptomic profiling identifies a high-stemness malignant-cell program associated with PI3K-AKT-mTOR signaling in lung adenocarcinoma"*.

## Study overview
We integrated single-cell RNA sequencing, spatial transcriptomics and TCGA-LUAD bulk data to define a high-stemness (HighStem) malignant-cell program in lung adenocarcinoma (LUAD).
- **Malignant cells:** 13,425 malignant epithelial cells, identified with inferCNV and scored with CytoTRACE, were split by quartiles into LowStem, DTStem and HighStem groups.
- **Characterisation:** the HighStem state was characterised with scPagwas tumour relevance scores, CellChat communication analysis and pathway activity.
- **Activity gene set:** hdWGCNA hub genes were filtered by Pearson correlation with CytoTRACE score within HighStem cells. This gave a 138-gene HighStem activity gene set, which was projected onto a spatial section and evaluated in TCGA-LUAD.
- **Six-gene panel:** six feature-selection methods retained six genes (TMSB10, S100A16, HLA-C, S100A14, CKLF, SRPX2).
- **Follow-up of the panel:** it was evaluated as an internal HighStem-versus-LowStem classifier, gene by gene in TCGA-LUAD, and for correlation with PI3K-AKT-mTOR axis genes.

## Data sources
| Data | Source | Expected local path |
|---|---|---|
| Single-cell RNA-seq of LUAD and normal lung | GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | `data/raw/GSE131907/`; processed objects in `data/processed/` |
| Spatial transcriptomics (P20_LUAD, Visium CytAssist) | GEO [GSE307534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307534) | `data/spatial/GSE307534/` |
| TCGA-LUAD expression and survival | NCI Genomic Data Commons / UCSC Xena ([TCGA-LUAD](https://portal.gdc.cancer.gov/projects/TCGA-LUAD)) | `data/external/tcga_luad/expression.tsv.gz`, `survival.tsv.gz`; Ensembl-to-symbol map `data/external/ensembl_to_symbol.csv` |

No data are distributed with this repository. All cell-level and patient-level data must be obtained from the sources above.

## Directory conventions
- **Working directory:** run every script from the repository root. Set `PROJECT_ROOT` to the repository root for the shell drivers, e.g. `export PROJECT_ROOT=$(pwd)`.
- **Paths:** scripts read inputs from `data/` and write outputs under `results/reproduce_paper/runs/<run name>/`. All paths are relative to the repository root.
- **Intermediate run directories:** several scripts read tables written by earlier steps. Their names are the run directories used in the original analysis and are kept unchanged so that scripts and outputs match, for example:
  - `paper_main_epiref_latest13k_originalflow_20260214_125353`: hdWGCNA network.
  - `final_malignant_v3_20261009/J2`: classifier cache.
  - `final_malignant_v3_20261009/J4/main_1200`: CellChat main analysis.
- **Script name suffixes:** files ending in `_copy` are minimally modified versions of the corresponding original script and are the versions that produced the published results. The original script is kept alongside where it is a direct parent.

## Run order and figure mapping
| Step | Folder | Main scripts | Produces |
|---|---|---|---|
| 1 | `01_atlas_and_stemness/` | `fig1_seurat_harmony.R` (Seurat SCTransform, Harmony, clustering at resolution 0.8, annotation); `00_export_for_seurat.py`; `01_c_cytotrace_manual.py` (CytoTRACE score, see below) | single-cell atlas; CytoTRACE scores |
| 2 | `02_malignant_calling_trs/` | `00c_epi_ref_main_step17.R` (inferCNV); `00c4_…`, `00c5_…` (malignant calling, reference mean + 3 SD); `00d1_…` (quartile stemness groups); `13_scpagwas_full.R`, `13c_trs_full_from_15k.R` (TRS); `prepare_latest13k_tables.py`; `00c8_…`, `00c8a_…` (original Fig. 2 panels) | malignant-cell metadata (13,425 cells; 3,357 / 6,711 / 3,357) |
| 3 | `03_cellchat/` | `run_J4_main_and_sens.sh` (calls `18a_…`, `18b_…` and the two sensitivity copies); `run_cellchat_tme_by_sample_planB_J4_copy.sh` + `18c_…` (sample-level tests); `partner_level_recompute.py`, `analyze_cellchat_runs.py`, `make_corrected_summary.py` | CellChat objects and tables (Fig. 3, Fig. S3, Fig. 8A) |
| 4 | `04_hdwgcna_feature_selection/` | `02b_run_hdwgcna_from_newmal_rds.R` (hdWGCNA, metacells k = 50, max_shared = 15, min_cells = 51, soft power 10); `03_ml_screening_final_useRaw_scores_copy.py` (Pearson filter r > 0.15, p < 0.05; six-method consensus); `41_native_ml_consensus_from_run.R` (native-R robustness check); `summarize_I1.py` | 334 hub genes; 138 activity genes; six-gene panel |
| 5 | `05_activity_signature/` | `export_fig5_activity_scores.py` (TCGA scores); `26_spatial_module_score.R` / `score_spatial_copy.R` (spatial AddModuleScore); `40_…` (tumour vs normal); `35_…` (KM, optimal cutpoint); `fig5D_three_ways.py` | Fig. 5 source tables |
| 6 | `06_classifier/` | `run_J2.sh` (repeated CV and full panel set with scikit-learn class labels); `make_fig6_reference_like_FORCE_HGB_copy.py` (HistGradientBoosting for fixed-split panels); `make_fig6_reference_like_PANELC_13664_copy.py`, `panelC_auc_by_gene.py`, `make_panelC_curves_P3.py` (single-gene ROC); `export_fig6_tables_six.py` | classifier cache and Fig. 6 tables |
| 7 | `07_tcga_genes/` | `38_make_fig7_SIXGENE_O5_copy.py` (data loading and statistics for the six genes); `supp_gene_stats_six.py` | Fig. 7 statistics |
| 8 | `08_pi3k_axis/` | `make_fig8_bioinfo_pi3k_allraw_SIXGENE_primary517_textfix_copy.R` (517 primary tumours, Spearman) | Fig. 8 source tables |
| 9 | `figures/` | one folder per figure (below); `figures/common/figlib.py` (180 mm layout helpers and checks); `figures/final_package/build_final.py` (PDF/TIF export and checks) | final figures |

| Figure | Script(s) in `figures/` | Panels |
|---|---|---|
| Fig. 1 | `Fig1/Fig1_export_data.R` → `Fig1/make_Fig1.py` | A UMAP clusters, B cell types, C marker dot plot |
| Fig. 2 | `Fig2/make_Fig2.py` | A CytoTRACE UMAP, B score distribution, C stemness groups, D CytoTRACE by group, E TRS by group, F CytoTRACE vs TRS |
| Fig. 3 | `Fig3/Fig3D_heatmap.R` → `Fig3/make_Fig3.py` | A outgoing vs incoming strength, B–C partner-level pair counts, D pathway activity |
| Fig. 4 | `Fig4/Fig4_panels.R` → `Fig4/make_Fig4.py` | A dendrogram, B kME, C module UMAP, D module dot plot |
| Fig. 5 | `Fig5/Fig5B_spatial.R` → `Fig5/make_Fig5.py` | A activity genes, B spatial projection, C TCGA tumour vs normal, D KM, E UpSet |
| Fig. 6 | `Fig6/make_Fig6.py` | A single-gene ROC, B repeated CV, C development/test ROC and PR, D confusion matrices, E SHAP |
| Fig. 7 | `Fig7/fig7lib.py`, `Fig7/make_Fig7.py` | A tumour vs normal, B paired samples, C ROC, D S100A16 KM |
| Fig. 8 | `Fig8/make_Fig8.py` | A pathway activity, B gene × PI3K-axis rho, C program score vs PI3K-axis genes |
| Fig. S1 | `FigS1/FigS1_density.R` → `FigS1/make_FigS1.py` | marker density maps |
| Fig. S2 | `FigS2/FigS2A_cnv_heatmap.R` → `FigS2/make_FigS2.py` | CNV heatmap, CNV score, TRS UMAP, stem/TRS overlay |
| Fig. S3 | `FigS3/FigS3AB_cellchat.R` → `FigS3/make_FigS3.py` | CellChat circles, role heatmaps, ligand-receptor bubbles, metabolic pathways |
| Fig. S4 | `FigS4/FigS4_panels.R` → `FigS4/make_FigS4.py` | soft-threshold diagnostics, module correlation |
| Fig. S5 | `FigS5/make_FigS5.py` | Boruta, decision tree, random forest, GBM, LASSO, ABESS proxy |
| Fig. S6 | `FigS6/FigS6_export_density_data.py` → `FigS6/make_FigS6.py` | density maps, model PR/ROC, decision curve, dependence plots |
| Fig. S7 | `FigS7/make_FigS7.py` | KM curves of the other five genes |
| Fig. S8 | `FigS8/make_FigS8.py` | six genes vs PDPK1/AKT1/MTOR |

## Environment
- **R 4.5.2** with the packages in `environment/R_packages.txt`. abess, Boruta, gbm and xgboost were installed in a project library `.Rlibs/`; the scripts that need them put it first on the library path.
- **Python 3.11.14** for the main analyses (`environment/requirements_main.txt`).
- **Python 3.12.12** for the Fig. 6 classifier (`environment/requirements_classifier.txt`). Python drivers read the interpreter from `PY_BIN`.
- **External tools:** poppler (`pdftotext`, `pdftohtml`) is used by the figure checks.
- **Random seeds:**
  - Seed 42 for feature selection, classifier splits and CellChat sampling (Boruta round *i* uses 42 + *i*).
  - Seed 20260125 for the Fig. 1 atlas.

## Software versions
| Software | Version | Used for |
|---|---|---|
| R | 4.5.2 | all R analyses |
| Seurat (SeuratObject) | 5.4.0 (5.3.0) | QC, normalisation, clustering; spatial module score |
| harmony | 1.2.4 | integration |
| infercnv | 1.26.0 | malignant-cell calling |
| scPagwas | 2.0.0 | tumour relevance score |
| CellChat | 1.6.1 | cell-cell communication |
| hdWGCNA / WGCNA | 0.4.09 / 1.73 | co-expression network |
| survival / survminer | 3.8-3 / 0.5.1 | survival analysis (Fig. 5D) |
| abess / Boruta | 0.4.11 / 9.0.0 | native-R robustness check only |
| ComplexHeatmap / ggplot2 | 2.26.0 / 4.0.2 | figures |
| Python | 3.11.14 / 3.12.12 | main analyses / classifier |
| scikit-learn | 1.8.0 | feature selection, classifiers |
| shap | 0.51.0 | SHAP values |
| numpy | 1.26.4 / 2.4.3 | general |
| pandas | 2.3.3 | general |
| scipy | 1.16.3 / 1.17.1 | statistical tests |
| scanpy / anndata | 1.11.5 / 0.12.9 (0.12.10) | h5ad handling |
| lifelines | 0.30.0 | survival analysis (Fig. 7) |
| matplotlib | 3.10.8 | figures |

More details are in `environment/software_versions.md`.

## Implementation notes
- **CytoTRACE:** no CytoTRACE package was used. `01_atlas_and_stemness/01_c_cytotrace_manual.py` implements the CytoTRACE gene-counts approach in Python: take the 200 genes most correlated with the number of detected genes, average their expression, smooth over a 30-nearest-neighbour graph, then rank-normalise. Scores were computed once for 33,697 epithelial cells; stemness quartiles were then taken within the 13,425 malignant cells.
- **Python feature selection:** Boruta is a simplified implementation (20 rounds of random forest with shadow features and a binomial test). "ABESS" is a proxy, recursive feature elimination with a linear SVM keeping 15 genes. LASSO is Gaussian `LassoCV`.
- **Native-R robustness check:** `41_native_ml_consensus_from_run.R` uses the abess and Boruta R packages.
- **Fig. 6 evaluation:** panel B is repeated stratified 5-fold cross-validation (10 repeats); panels C–E use one fixed stratified split (test size 0.2, seed 42).

## License
MIT (see `LICENSE`).
