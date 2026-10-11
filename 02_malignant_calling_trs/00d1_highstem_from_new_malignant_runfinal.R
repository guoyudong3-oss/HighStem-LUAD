#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(Seurat)
  library(future)
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

must_num <- function(x, name) {
  xn <- suppressWarnings(as.numeric(x))
  if (is.na(xn)) stop("Invalid numeric for ", name, ": ", x)
  xn
}

must_int <- function(x, name) {
  xi <- suppressWarnings(as.integer(x))
  if (is.na(xi) || xi <= 0) stop("Invalid integer for ", name, ": ", x)
  xi
}

args <- commandArgs(trailingOnly = TRUE)
run_dir <- get_arg(args, "--run-dir", NULL)
method <- get_arg(args, "--method", "ref_mean_3sd")
q_low <- must_num(get_arg(args, "--q-low", "0.25"), "--q-low")
q_high <- must_num(get_arg(args, "--q-high", "0.75"), "--q-high")
min_pct <- must_num(get_arg(args, "--min-pct", "0.10"), "--min-pct")
logfc_th <- must_num(get_arg(args, "--logfc-threshold", "0.25"), "--logfc-threshold")
workers <- must_int(get_arg(args, "--workers", "6"), "--workers")

cyto_meta_file <- get_arg(args, "--cytotrace-meta", "results/reproduce_paper/r_data/metadata.csv")
seurat_rds <- get_arg(args, "--seurat-rds", "data/processed/step1_seurat_harmony.rds")

if (is.null(run_dir) || !nzchar(run_dir)) {
  stop("Usage: Rscript 00d1_highstem_from_new_malignant_runfinal.R --run-dir <path> [--method ref_mean_3sd] [--q-low 0.25] [--q-high 0.75] [--workers 6]")
}

run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
tables_dir <- file.path(run_dir, "tables")
fig_dir <- file.path(run_dir, "figures")
rds_dir <- file.path(run_dir, "rds")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

malignant_meta_file <- file.path(tables_dir, paste0("cnv_paper_exact_cell_metadata_runfinal_", method, ".csv"))
if (!file.exists(malignant_meta_file)) stop("Missing: ", malignant_meta_file)
if (!file.exists(cyto_meta_file)) stop("Missing: ", cyto_meta_file)
if (!file.exists(seurat_rds)) stop("Missing: ", seurat_rds)

save_plot_both <- function(p, base, w, h, dpi = 320) {
  ggsave(file.path(fig_dir, paste0(base, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(base, ".png")), p, width = w, height = h, dpi = dpi)
}

message("Loading malignant calls: ", malignant_meta_file)
meta <- fread(malignant_meta_file)
meta_m <- meta[infercnv_group == "Observation" & malignant_status == "malignant"]
if (nrow(meta_m) < 100) stop("Too few malignant observation cells: ", nrow(meta_m))
message("Malignant observation cells: ", nrow(meta_m))

message("Loading CytoTRACE metadata: ", cyto_meta_file)
cyto <- fread(cyto_meta_file, select = c("cell_id", "cytotrace_score"))

df <- merge(meta_m, cyto, by = "cell_id", all.x = TRUE)
missing_ct <- sum(is.na(df$cytotrace_score))
message("Malignant cells missing cytotrace_score: ", missing_ct)

df <- df[!is.na(cytotrace_score)]
if (nrow(df) < 200) stop("Too few malignant cells with cytotrace_score available: ", nrow(df))

qs <- quantile(df$cytotrace_score, probs = c(q_low, q_high), na.rm = TRUE)
cut_low <- as.numeric(qs[[1]])
cut_high <- as.numeric(qs[[2]])

df[, stem_group := ifelse(cytotrace_score <= cut_low, "LowStem",
                          ifelse(cytotrace_score >= cut_high, "HighStem", "DTStem"))]
df[, stem_group := factor(stem_group, levels = c("LowStem", "DTStem", "HighStem"))]

out_meta <- file.path(tables_dir, paste0("stemness_cytotrace_malignant_runfinal_", method, ".csv"))
fwrite(df[, .(cell_id, Sample_Origin, infercnv_group, cnv_score_runfinal, malignant_status, cytotrace_score, stem_group)], out_meta)
message("Saved: ", out_meta)

summary_file <- file.path(tables_dir, paste0("stemness_cytotrace_summary_malignant_runfinal_", method, ".txt"))
counts <- table(df$stem_group)
txt <- c(
  paste0("Run dir: ", run_dir),
  paste0("Malignant calling: run.final + ", method),
  paste0("CytoTRACE source: ", normalizePath(cyto_meta_file, winslash = "/")),
  paste0("Malignant obs cells (total): ", nrow(meta_m)),
  paste0("Malignant obs cells (with cytotrace_score): ", nrow(df)),
  paste0("Missing cytotrace_score: ", missing_ct),
  paste0("Quantile cutoffs: q_low=", q_low, " (", format(cut_low, digits = 6), "), q_high=", q_high, " (", format(cut_high, digits = 6), ")"),
  "Stem group counts:",
  paste(names(counts), as.integer(counts), sep = ": ")
)
writeLines(txt, summary_file)
message("Saved: ", summary_file)

message("Plotting CytoTRACE distribution + stem group violin ...")
theme_fig <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))
stem_cols <- c("LowStem" = "#56B1A6", "DTStem" = "#C7C7C7", "HighStem" = "#C92D45")

p_hist <- ggplot(df, aes(x = cytotrace_score)) +
  geom_histogram(bins = 60, fill = "grey85", color = "grey35", linewidth = 0.2) +
  geom_vline(xintercept = c(cut_low, cut_high), linetype = "dashed", color = "#C92D45") +
  labs(title = "CytoTRACE score distribution (malignant)", x = "CytoTRACE score", y = "Cells") +
  theme_fig
save_plot_both(p_hist, paste0("Stem_CytoTRACE_hist_malignant_", method), 6.2, 3.6)

p_vln <- ggplot(df, aes(x = stem_group, y = cytotrace_score, fill = stem_group)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.18, outlier.size = 0.2, alpha = 0.9) +
  scale_fill_manual(values = stem_cols) +
  labs(title = "CytoTRACE by stem group (malignant)", x = NULL, y = "CytoTRACE score") +
  theme_fig +
  theme(legend.position = "none")
save_plot_both(p_vln, paste0("Stem_CytoTRACE_violin_malignant_", method), 4.6, 3.9)

message("Loading Seurat: ", seurat_rds)
seu <- readRDS(seurat_rds)

cells_use <- intersect(colnames(seu), df$cell_id)
if (length(cells_use) < 200) stop("Too few cells found in Seurat object after intersection: ", length(cells_use))
df2 <- df[match(cells_use, df$cell_id)]

seu_m <- subset(seu, cells = cells_use)
seu_m$cytotrace_score <- df2$cytotrace_score
seu_m$stem_group <- df2$stem_group

out_seu <- file.path(rds_dir, paste0("seurat_malignant_with_stemness_runfinal_", method, ".rds"))
saveRDS(seu_m, out_seu)
message("Saved: ", out_seu)

message("UMAP plots (malignant only, using existing Seurat UMAP) ...")
um <- Embeddings(seu_m, "umap")
plot_df <- data.frame(
  cell_id = rownames(um),
  UMAP_1 = um[, 1],
  UMAP_2 = um[, 2],
  cytotrace_score = seu_m$cytotrace_score,
  stem_group = seu_m$stem_group,
  stringsAsFactors = FALSE
)

p_um_ct <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = cytotrace_score)) +
  geom_point(size = 0.12, alpha = 0.9) +
  scale_color_gradientn(colours = c("#2c7bb6", "#ffffbf", "#d7191c")) +
  labs(title = "CytoTRACE (malignant)", x = "UMAP_1", y = "UMAP_2", color = "Score") +
  theme_fig
save_plot_both(p_um_ct, paste0("Stem_CytoTRACE_umap_malignant_", method), 5.4, 4.4)

p_um_grp <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = stem_group)) +
  geom_point(size = 0.12, alpha = 0.9) +
  scale_color_manual(values = stem_cols) +
  labs(title = "Stem groups (malignant)", x = "UMAP_1", y = "UMAP_2", color = "Group") +
  theme_fig
save_plot_both(p_um_grp, paste0("Stem_group_umap_malignant_", method), 5.4, 4.4)

message("DEG: HighStem vs LowStem (wilcox) ...")
seu_de <- subset(seu_m, subset = stem_group %in% c("HighStem", "LowStem"))
Idents(seu_de) <- "stem_group"

future::plan("multisession", workers = workers)
options(future.globals.maxSize = 10 * 1024^3)

deg <- FindMarkers(
  seu_de,
  ident.1 = "HighStem",
  ident.2 = "LowStem",
  test.use = "wilcox",
  min.pct = min_pct,
  logfc.threshold = logfc_th
)

deg$gene <- rownames(deg)
deg <- deg[, c("gene", setdiff(colnames(deg), "gene"))]
deg_out <- file.path(tables_dir, paste0("deg_highstem_vs_lowstem_malignant_runfinal_", method, ".csv"))
fwrite(as.data.table(deg), deg_out)
message("Saved: ", deg_out)

message("Done.")

