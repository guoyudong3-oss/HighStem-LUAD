#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

get_arg <- function(args, key, default = NULL) {
  key_eq <- paste0("^", key, "=")
  idx <- which(args == key)
  if (length(idx) > 0) {
    if (idx[length(idx)] >= length(args)) stop("Missing value for ", key)
    return(args[idx[length(idx)] + 1])
  }
  kv <- grep(key_eq, args, value = TRUE)
  if (length(kv) > 0) return(sub(key_eq, "", kv[length(kv)]))
  default
}

must_int <- function(x, name) {
  xi <- suppressWarnings(as.integer(x))
  if (is.na(xi) || xi <= 0) stop("Invalid integer for ", name, ": ", x)
  xi
}

args <- commandArgs(trailingOnly = TRUE)
run_dir <- get_arg(args, "--run-dir", NULL)
method <- get_arg(args, "--method", "ref_mean_3sd")
k_clusters <- must_int(get_arg(args, "--k", "4"), "--k")
genes_per_bin <- must_int(get_arg(args, "--genes-per-bin", "25"), "--genes-per-bin")
group_mode <- get_arg(args, "--group-mode", "malignant_status") # malignant_status|none

if (is.null(run_dir) || !nzchar(run_dir)) {
  stop("Usage: Rscript 00c8a_make_fig2A_like_paper_runfinal.R --run-dir <path> [--method ref_mean_3sd] [--k 4] [--genes-per-bin 25] [--group-mode malignant_status|none]")
}

run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
tables_dir <- file.path(run_dir, "tables")
fig_dir <- "results/reproduce_paper/runs/final_malignant_v3_20261009/W/FigS2/A/out"
infer_dir <- file.path(run_dir, "infercnv_output")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


infer_rds <- file.path(infer_dir, "run.final.infercnv_obj")
if (!file.exists(infer_rds)) stop("Missing: ", infer_rds)

meta_malignant_file <- file.path(tables_dir, paste0("cnv_paper_exact_cell_metadata_runfinal_", method, ".csv"))
if (group_mode == "malignant_status" && !file.exists(meta_malignant_file)) {
  stop("Missing: ", meta_malignant_file, "\nExpected this file to provide malignant_status for group annotation.")
}

message("Loading inferCNV run.final ...")
infer_obj <- readRDS(infer_rds)
expr <- infer_obj@expr.data
go <- infer_obj@gene_order

if (!is.matrix(expr)) stop("infer_obj@expr.data is not a matrix.")
if (!is.data.frame(go)) stop("infer_obj@gene_order is not a data.frame.")
if (!all(rownames(go) == rownames(expr))) {
  go <- go[match(rownames(expr), rownames(go)), , drop = FALSE]
}

obs_idx <- unlist(infer_obj@observation_grouped_cell_indices, use.names = FALSE)
obs_ids <- colnames(expr)[obs_idx]
if (length(obs_ids) == 0) stop("No observation cells found in inferCNV object.")
message("Observation cells: ", length(obs_ids))

expr_obs <- expr[, obs_ids, drop = FALSE] # genes x obs_cells

message("Gene-binning CNV profiles (", genes_per_bin, " genes/bin) ...")
chr <- as.character(go$chr)
bin_num <- integer(length(chr))
for (c in unique(chr)) {
  idx <- which(chr == c)
  bin_num[idx] <- ceiling(seq_along(idx) / genes_per_bin)
}
bin_id <- paste0(chr, "_", sprintf("%04d", bin_num))
bin_levels <- unique(bin_id) # keep genomic order
bin_factor <- factor(bin_id, levels = bin_levels)

bin_sum <- rowsum(expr_obs, group = bin_factor, reorder = FALSE) # bins x cells
bin_counts <- as.vector(table(bin_factor))
names(bin_counts) <- levels(bin_factor)
bin_mean <- bin_sum / bin_counts[match(rownames(bin_sum), names(bin_counts))]
rm(bin_sum)
gc()

rm(expr_obs)
gc()

message("Using existing k-cluster assignment (same table as Fig. 2B) ...")
kc <- read.csv(file.path(tables_dir, paste0("fig2_k_cluster_runfinal_k", k_clusters, "_g", genes_per_bin, "_", method, ".csv")), stringsAsFactors = FALSE)
obs_ids <- obs_ids[obs_ids %in% kc$cell_id]
bin_mean <- bin_mean[, obs_ids, drop = FALSE]
df_k <- data.frame(cell_id = obs_ids, cnv_score_runfinal = kc$cnv_score_runfinal[match(obs_ids, kc$cell_id)], stringsAsFactors = FALSE)
df_k$k_cluster <- factor(as.character(kc$k_cluster[match(obs_ids, kc$cell_id)]), levels = as.character(seq_len(k_clusters)))
rm(feat)
message("Panel A (cell-level, paper-like): CNV heatmap ...")
ord <- order(df_k$k_cluster, df_k$cnv_score_runfinal, decreasing = FALSE)
ord_ids <- df_k$cell_id[ord]

mat <- t(bin_mean[, ord_ids, drop = FALSE]) # cells x bins
rm(bin_mean)
gc()

bin_cols <- colnames(mat)
bin_chr <- sub("(_.*)$", "", bin_cols)
chr_levels <- unique(bin_chr)
chr_fac <- factor(gsub("^chr", "", bin_chr), levels = gsub("^chr", "", chr_levels))

vmin <- as.numeric(quantile(mat, 0.01, na.rm = TRUE))
vmax <- as.numeric(quantile(mat, 0.99, na.rm = TRUE))
col_fun <- circlize::colorRamp2(c(vmin, 1, vmax), c("#2166ac", "#f7f7f7", "#b2182b"))

k_pal <- c("1" = "#8DD3C7", "2" = "#FFFFB3", "3" = "#B3DE69", "4" = "#FB8072")
clusters <- levels(df_k$k_cluster)

group_lab <- NULL
if (group_mode == "malignant_status") {
  meta <- read.csv(meta_malignant_file, stringsAsFactors = FALSE)
  ms <- meta$malignant_status[match(rownames(mat), meta$cell_id)]
  if (any(is.na(ms))) {
    n_na <- sum(is.na(ms))
    warning("Missing malignant_status for ", n_na, " cells; treating as Normal.")
    ms[is.na(ms)] <- "normal_epithelial"
  }
  group_lab <- ifelse(ms == "malignant", "Tumor", "Normal")
  group_lab <- factor(group_lab, levels = c("Normal", "Tumor"))
}

ra <- rowAnnotation(
  k_cluster = df_k$k_cluster[match(rownames(mat), df_k$cell_id)],
  col = list(
    k_cluster = k_pal[clusters],
    group = c("Normal" = "#30A3A1", "Tumor" = "#D73027")
  ),
  annotation_name_gp = gpar(fontsize = 8, fontfamily = "Arial"),
  width = unit(6, "mm")
)

if (!is.null(group_lab)) {
  ra <- rowAnnotation(
    k_cluster = df_k$k_cluster[match(rownames(mat), df_k$cell_id)],
    group = group_lab,
    col = list(
      k_cluster = k_pal[clusters],
      group = c("Normal" = "#30A3A1", "Tumor" = "#D73027")
    ),
    annotation_name_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    width = unit(7, "mm")
  )
}

ht <- Heatmap(
  mat,
  name = "CNV",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_quality = 3,
  column_split = chr_fac,
  column_gap = unit(0.2, "mm"),
  column_title_gp = gpar(fontsize = 8, fontfamily = "Arial"), column_title_rot = 90,
  row_split = df_k$k_cluster[match(rownames(mat), df_k$cell_id)],
  row_gap = unit(0.6, "mm"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 8, fontfamily = "Arial")
  )
)

# Use cairo_pdf to ensure embedded TrueType fonts (AI-friendly).
cairo_pdf(file.path(fig_dir, "FigS2A_cnv_heatmap_v2.pdf"), width = 7.0866, height = 2.75, family = "Arial")
draw(ra + ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
message("Done")
