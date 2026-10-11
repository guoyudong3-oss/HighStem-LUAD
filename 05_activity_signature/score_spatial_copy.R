# Scoring steps copied from 26_spatial_module_score.R: Read10X -> CreateSeuratObject -> NormalizeData -> AddModuleScore (default parameters)
# Difference: gene set and output path are given as arguments.
suppressPackageStartupMessages({library(Seurat)})
args <- commandArgs(trailingOnly = TRUE)
sample_path <- args[1]; gene_file <- args[2]; out_csv <- args[3]
activity_genes <- read.csv(gene_file)$gene
sample_name <- basename(sample_path)
counts_raw <- Read10X(file.path(sample_path, "filtered_feature_bc_matrix"))
counts <- if (is.list(counts_raw)) { if ("Gene Expression" %in% names(counts_raw)) counts_raw$`Gene Expression` else counts_raw[[1]] } else counts_raw
pos_path <- file.path(sample_path, "spatial", "tissue_positions.csv")
if (!file.exists(pos_path)) pos_path <- file.path(sample_path, "spatial", "tissue_positions_list.csv")
coords <- read.csv(pos_path, header = FALSE)
if (!grepl("-1", coords[1, 1])) coords <- read.csv(pos_path, header = TRUE)
colnames(coords)[1:6] <- c("barcode", "tissue", "row", "col", "imagerow", "imagecol")
rownames(coords) <- coords$barcode
common <- intersect(colnames(counts), rownames(coords))
st_obj <- CreateSeuratObject(counts = counts[, common], assay = "Spatial")
st_obj <- NormalizeData(st_obj, verbose = FALSE)
valid_genes <- intersect(activity_genes, rownames(st_obj))
st_obj <- AddModuleScore(st_obj, features = list(valid_genes), name = "HighStem_Activity_Score")
# Independent epithelial/tumour marker score (signature genes excluded), used only to define regions
epi <- intersect(c("EPCAM","KRT7","KRT8","KRT18","KRT19","NAPSA","NKX2-1","SFTPB"), rownames(st_obj))
epi <- setdiff(epi, activity_genes)
st_obj <- AddModuleScore(st_obj, features = list(epi), name = "EpiMarker_Score")
res_df <- data.frame(Barcode = colnames(st_obj), Score = st_obj$HighStem_Activity_Score1,
  EpiMarker_Score = st_obj$EpiMarker_Score1,
  X_pxl = coords[colnames(st_obj), "imagecol"], Y_pxl = coords[colnames(st_obj), "imagerow"],
  In_Tissue = coords[colnames(st_obj), "tissue"])
write.csv(res_df, out_csv, row.names = FALSE)
cat(sprintf("SAMPLE=%s spots=%d genes_in_set=%d genes_valid=%d epi_markers=%s\n", sample_name, nrow(res_df), length(activity_genes), length(valid_genes), paste(epi, collapse=",")))
cat("MISSING:", paste(setdiff(activity_genes, rownames(st_obj)), collapse=","), "\n")
