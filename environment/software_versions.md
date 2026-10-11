# Software versions

Read from the environments used for the analyses (R: Homebrew R with the project library `.Rlibs/` first on the library path; Python: the project virtual environments).

## R 4.5.2
| Package | Version | Used for |
|---|---|---|
| Seurat (SeuratObject) | 5.4.0 (5.3.0) | QC, SCTransform, dimensionality reduction, clustering (Fig. 1); spatial AddModuleScore (Fig. 5B) |
| harmony | 1.2.4 | sample integration (Fig. 1) |
| infercnv | 1.26.0 | copy-number inference and malignant-cell calling (Fig. 2, Fig. S2) |
| scPagwas | 2.0.0 | tumour relevance score (TRS; gene weights estimated on 15,000 randomly sampled cells, TRS computed for all cells) |
| CellChat | 1.6.1 | cell-cell communication (Fig. 3, Fig. S3) |
| hdWGCNA | 0.4.09 | co-expression network (Fig. 4, Fig. S4) |
| WGCNA | 1.73 | soft-threshold selection, TOM (Fig. 4, Fig. S4) |
| survival | 3.8-3 | Cox models, log-rank tests (Fig. 5D) |
| survminer (maxstat) | 0.5.1 (0.7-26) | optimal cutpoint and KM curves (Fig. 5D) |
| abess | 0.4.11 | native-R robustness check of feature selection only |
| Boruta | 9.0.0 | native-R robustness check of feature selection only |
| ComplexHeatmap | 2.26.0 | heatmaps (Fig. 3D, Fig. S2A, Fig. S3B) |
| ggplot2 | 4.0.2 | R figure panels |

## Python
| Package | Main environment (Python 3.11.14) | Classifier environment (Python 3.12.12) | Used for |
|---|---|---|---|
| scikit-learn | 1.8.0 | 1.8.0 | six-method feature selection; classifiers and cross-validation (Fig. 6) |
| shap | not used for results | 0.51.0 | SHAP values (Fig. 6) |
| numpy | 1.26.4 | 2.4.3 | general |
| pandas | 2.3.3 | 2.3.3 | general |
| scipy | 1.16.3 | 1.17.1 | statistical tests |
| scanpy | 1.11.5 | not installed | reading and handling h5ad files |
| anndata | 0.12.9 | 0.12.10 | reading h5ad files |
| lifelines | 0.30.0 | not installed | KM and Cox models (Fig. 7) |
| matplotlib | 3.10.8 | 3.10.8 | figures |

## CytoTRACE
CytoTRACE scores were computed with `01_atlas_and_stemness/01_c_cytotrace_manual.py`, a Python implementation of the CytoTRACE gene-counts approach (no CytoTRACE package was used): the 200 genes most positively correlated with the number of detected genes per cell are averaged, smoothed over a 30-nearest-neighbour graph and converted to rank percentiles. This step was run on 2026-01-26 in an earlier Python environment whose package versions were not recorded.
