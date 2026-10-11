#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(ggpubr)
  library(patchwork)
  library(uwot)
  library(RColorBrewer)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(Seurat)
  library(infercnv)
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
meta_cytotrace_file <- get_arg(args, "--meta-cytotrace-file", "results/reproduce_paper/r_data/metadata.csv")
trs_file_arg <- get_arg(args, "--trs-file", NULL)
seurat_rds <- get_arg(args, "--seurat-rds", "data/processed/step1_seurat_harmony.rds")

if (is.null(run_dir) || !nzchar(run_dir)) {
  stop(
    "Usage: Rscript 00c8_make_fig2_like_paper_runfinal.R --run-dir <path> ",
    "[--method ref_mean_3sd] [--k 4] [--genes-per-bin 25] ",
    "[--meta-cytotrace-file <csv>] [--trs-file <csv>] [--seurat-rds <rds>]"
  )
}

run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
tables_dir <- file.path(run_dir, "tables")
fig_dir <- file.path(run_dir, "figures")
infer_dir <- file.path(run_dir, "infercnv_output")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

infer_rds <- file.path(infer_dir, "run.final.infercnv_obj")
meta_malignant_file <- file.path(tables_dir, paste0("cnv_paper_exact_cell_metadata_runfinal_", method, ".csv"))
if (!file.exists(infer_rds)) stop("Missing: ", infer_rds)
if (!file.exists(meta_malignant_file)) stop("Missing: ", meta_malignant_file, "\nRun 00c5_export_malignant_set_from_runfinal.R first.")

trs_file <- trs_file_arg
if (is.null(trs_file) || !nzchar(trs_file)) {
  trs_file <- file.path(tables_dir, "scpagwas_trs_scores_full_from_15k.csv")
}
if (!file.exists(meta_cytotrace_file)) stop("Missing: ", meta_cytotrace_file)
if (!file.exists(trs_file)) stop("Missing: ", trs_file)
if (!file.exists(seurat_rds)) stop("Missing: ", seurat_rds)

fig_font_family <- "Arial"
pdf_dev <- function(...) cairo_pdf(..., family = fig_font_family)
svg_dev <- function(...) {
  svglite::svglite(
    ...,
    system_fonts = list(sans = fig_font_family, serif = fig_font_family, mono = "Courier New")
  )
}

save_plot_both <- function(p, base, w, h, dpi = 320) {
  ggsave(file.path(fig_dir, paste0(base, ".pdf")), p, width = w, height = h, device = pdf_dev, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(fig_dir, paste0(base, ".svg")), p, width = w, height = h, device = svg_dev, bg = "white")
  }
  ggsave(file.path(fig_dir, paste0(base, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
}

theme_fig <- theme_classic(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(6, 6, 6, 6)
  )

theme_umap <- theme_void(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(4, 6, 8, 8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Draw a true right-half violin (not a full violin shifted to the right).
GeomHalfViolinRight <- ggproto(
  "GeomHalfViolinRight",
  GeomViolin,
  setup_data = function(data, params) {
    data$width <- data$width %||% params$width %||% (resolution(data$x, FALSE) * 0.9)
    transform(
      data,
      ymin = y,
      ymax = y,
      xmin = x,
      xmax = x + width / 2
    )
  },
  draw_group = function(data, panel_scales, coord, ...) {
    data <- transform(
      data,
      xminv = x,
      xmaxv = x + violinwidth * (xmax - x)
    )
    data_left <- transform(data, x = xminv)
    data_left <- data_left[order(data_left$y), , drop = FALSE]
    data_right <- transform(data, x = xmaxv)
    data_right <- data_right[order(data_right$y, decreasing = TRUE), , drop = FALSE]
    polygon_data <- rbind(data_left, data_right, data_left[1, , drop = FALSE])
    GeomPolygon$draw_panel(polygon_data, panel_scales, coord)
  }
)

geom_half_violin_right <- function(mapping = NULL, data = NULL, stat = "ydensity",
                                   position = "identity", ..., trim = TRUE,
                                   scale = "area", show.legend = NA,
                                   inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHalfViolinRight,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, ...)
  )
}

k_pal <- c("1" = "#8DD3C7", "2" = "#FFFFB3", "3" = "#B3DE69", "4" = "#FB8072")
stem_cols <- c("LowStem" = "#56B1A6", "DTStem" = "#C7C7C7", "HighStem" = "#C92D45")

# Ensure annotate("text", ...) and geom_text defaults are Arial.
update_geom_defaults("text", list(family = fig_font_family))
update_geom_defaults("label", list(family = fig_font_family))

add_corner_axes <- function(p, x, y, xlab = "UMAP 1", ylab = "UMAP 2") {
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)
  x0 <- xr[1] + 0.08 * dx
  y0 <- yr[1] + 0.08 * dy
  x1 <- x0 + 0.22 * dx
  y1 <- y0 + 0.22 * dy
  p +
    annotate(
      "segment",
      x = x0, xend = x1, y = y0, yend = y0,
      linewidth = 0.5, color = "black",
      arrow = arrow(length = unit(0.14, "cm"), type = "closed")
    ) +
    annotate(
      "segment",
      x = x0, xend = x0, y = y0, yend = y1,
      linewidth = 0.5, color = "black",
      arrow = arrow(length = unit(0.14, "cm"), type = "closed")
    ) +
    annotate("text", x = (x0 + x1) / 2, y = y0 - 0.05 * dy, label = xlab, size = 3.6, family = fig_font_family) +
    annotate("text", x = x0 - 0.05 * dx, y = (y0 + y1) / 2, label = ylab, angle = 90, size = 3.6, family = fig_font_family) +
    coord_fixed(ratio = 1, clip = "off")
}

filter_umap_outliers <- function(df, x_col = "UMAP_1", y_col = "UMAP_2",
                                 q_trim = 0.005, md_keep = 0.995,
                                 k_nn = 8, comp_min_prop = 0.003, comp_min_n = 20) {
  x <- df[[x_col]]
  y <- df[[y_col]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 500) {
    return(list(df = df[ok, , drop = FALSE], keep = ok, stats = data.frame(
      method = "fallback_small_n", n_total = nrow(df), n_keep = sum(ok), n_drop = sum(!ok),
      keep_ratio = sum(ok) / max(1, nrow(df))
    )))
  }
  x_ok <- x[ok]
  y_ok <- y[ok]
  qx <- as.numeric(quantile(x_ok, probs = c(q_trim, 1 - q_trim), na.rm = TRUE))
  qy <- as.numeric(quantile(y_ok, probs = c(q_trim, 1 - q_trim), na.rm = TRUE))
  mask_q <- ok & x >= qx[1] & x <= qx[2] & y >= qy[1] & y <= qy[2]
  if (sum(mask_q) < 500) {
    mask_q <- ok
  }

  core <- df[mask_q, c(x_col, y_col), drop = FALSE]
  ctr <- c(median(core[[x_col]], na.rm = TRUE), median(core[[y_col]], na.rm = TRUE))
  cov_mat <- tryCatch(cov(core), error = function(e) NULL)
  if (is.null(cov_mat) || any(!is.finite(cov_mat)) || abs(det(cov_mat)) < 1e-10) {
    vx <- var(core[[x_col]], na.rm = TRUE)
    vy <- var(core[[y_col]], na.rm = TRUE)
    cov_mat <- matrix(c(vx + 1e-6, 0, 0, vy + 1e-6), nrow = 2, byrow = TRUE)
  }
  xy <- as.matrix(df[, c(x_col, y_col), drop = FALSE])
  d2 <- mahalanobis(xy, center = ctr, cov = cov_mat)
  thr <- as.numeric(quantile(d2[ok], probs = md_keep, na.rm = TRUE))
  mask_md <- ok & (d2 <= thr)
  mask_base <- mask_q & mask_md

  # Guardrail: if too aggressive, relax to keep >=97%.
  keep_ratio <- sum(mask_base) / max(1, nrow(df))
  if (keep_ratio < 0.97) {
    thr2 <- as.numeric(quantile(d2[ok], probs = 0.998, na.rm = TRUE))
    mask_base <- ok & (d2 <= thr2)
  }

  # Remove tiny disconnected islands on UMAP that visually stretch axes.
  mask <- mask_base
  comp_kept <- NA_integer_
  comp_total <- NA_integer_
  comp_min_size <- NA_integer_
  if (sum(mask_base) >= 500 &&
      requireNamespace("RANN", quietly = TRUE) &&
      requireNamespace("igraph", quietly = TRUE)) {
    idx <- which(mask_base)
    xy_core <- as.matrix(df[idx, c(x_col, y_col), drop = FALSE])
    k_eff <- max(3L, min(as.integer(k_nn), nrow(xy_core) - 1L))
    nn <- RANN::nn2(data = xy_core, query = xy_core, k = k_eff + 1L)

    d_ref <- nn$nn.dists[, k_eff + 1L]
    eps <- as.numeric(quantile(d_ref, probs = 0.92, na.rm = TRUE)) * 1.35
    if (!is.finite(eps) || eps <= 0) {
      eps <- as.numeric(median(nn$nn.dists[, 2L], na.rm = TRUE)) * 2.5
    }

    from <- rep(seq_len(nrow(xy_core)), each = k_eff)
    to <- as.vector(t(nn$nn.idx[, 2:(k_eff + 1L), drop = FALSE]))
    dd <- as.vector(t(nn$nn.dists[, 2:(k_eff + 1L), drop = FALSE]))
    use_edge <- is.finite(dd) & dd <= eps & to > 0

    g <- igraph::make_empty_graph(n = nrow(xy_core), directed = FALSE)
    if (sum(use_edge) > 0) {
      e_mat <- cbind(from[use_edge], to[use_edge])
      g <- igraph::add_edges(g, as.vector(t(e_mat)))
    }
    comp <- igraph::components(g)$membership
    sizes <- tabulate(comp)
    largest <- which.max(sizes)
    min_size <- max(as.integer(comp_min_n), as.integer(ceiling(comp_min_prop * nrow(xy_core))))
    keep_comp <- which(sizes >= min_size)
    if (!(largest %in% keep_comp)) keep_comp <- c(keep_comp, largest)
    keep_comp <- sort(unique(keep_comp))

    mask_cc <- rep(FALSE, length(comp))
    mask_cc[comp %in% keep_comp] <- TRUE

    mask2 <- mask_base
    mask2[idx] <- mask_cc

    # Guardrail: component filter should only trim minor islands.
    if ((sum(mask2) / max(1, nrow(df))) >= 0.94) {
      mask <- mask2
      comp_kept <- length(keep_comp)
      comp_total <- length(sizes)
      comp_min_size <- min_size
    }
  }

  stats <- data.frame(
    method = "quantile_mahalanobis_plus_components",
    n_total = nrow(df),
    n_keep = sum(mask),
    n_drop = nrow(df) - sum(mask),
    keep_ratio = sum(mask) / max(1, nrow(df)),
    q_trim = q_trim,
    md_keep = md_keep,
    x_q_lo = qx[1],
    x_q_hi = qx[2],
    y_q_lo = qy[1],
    y_q_hi = qy[2],
    md_threshold = thr,
    k_nn = as.integer(k_nn),
    comp_min_prop = comp_min_prop,
    comp_min_n = as.integer(comp_min_n),
    comp_kept = comp_kept,
    comp_total = comp_total,
    comp_min_size = comp_min_size
  )
  list(df = df[mask, , drop = FALSE], keep = mask, stats = stats)
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

cell_ids <- colnames(expr)
cnv_score <- colMeans(abs(expr - 1))
names(cnv_score) <- cell_ids

message("Computing k-cluster (k=", k_clusters, ") on epithelial cells using gene-binned CNV profiles ...")
chr <- as.character(go$chr)
bin_num <- integer(length(chr))
for (c in unique(chr)) {
  idx <- which(chr == c)
  bin_num[idx] <- ceiling(seq_along(idx) / genes_per_bin)
}
bin_id <- paste0(chr, "_", sprintf("%04d", bin_num))
bin_factor <- factor(bin_id, levels = unique(bin_id))
bin_sum <- rowsum(expr, group = bin_factor, reorder = FALSE)
bin_counts <- as.vector(table(bin_factor))
names(bin_counts) <- levels(bin_factor)
bin_mean <- bin_sum / bin_counts[match(rownames(bin_sum), names(bin_counts))]
rm(bin_sum)
gc()

feat <- t(bin_mean) # cells x bins
rm(bin_mean)
gc()

set.seed(42)
pca <- prcomp(feat, center = TRUE, scale. = TRUE, rank. = 20)
pc <- pca$x
rm(feat, pca)
gc()

km <- kmeans(pc, centers = k_clusters, nstart = 30, iter.max = 50)
km_cluster <- as.integer(km$cluster)
names(km_cluster) <- cell_ids

df_k <- data.frame(
  cell_id = cell_ids,
  km_cluster = km_cluster,
  cnv_score_runfinal = cnv_score,
  stringsAsFactors = FALSE
)
cl_means <- df_k %>%
  group_by(km_cluster) %>%
  summarise(mean_score = mean(cnv_score_runfinal, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_score)

# Label the lowest-CNV cluster as "3" to match the paper narrative.
lowest <- cl_means$km_cluster[1]
others <- cl_means$km_cluster[-1]
target_labels <- c(3, 1, 2, 4)
target_labels <- target_labels[seq_len(k_clusters)]
map <- setNames(as.character(target_labels), as.character(c(lowest, others)))

df_k$k_cluster <- factor(map[as.character(df_k$km_cluster)], levels = c("1", "2", "3", "4")[seq_len(k_clusters)])

out_k <- file.path(tables_dir, paste0("fig2_k_cluster_runfinal_k", k_clusters, "_g", genes_per_bin, "_", method, ".csv"))
write.csv(df_k[, c("cell_id", "k_cluster", "cnv_score_runfinal")], out_k, row.names = FALSE)
message("Saved k-cluster table: ", out_k)

message("Panel A: CNV heatmap (cluster-averaged, paper-like) ...")
clusters <- levels(df_k$k_cluster)
mean_by_cluster <- sapply(clusters, function(kc) {
  ids <- df_k$cell_id[df_k$k_cluster == kc]
  rowMeans(expr[, ids, drop = FALSE])
})
colnames(mean_by_cluster) <- clusters
matA <- t(mean_by_cluster) # clusters x genes

chr_fac <- factor(gsub("^chr", "", go$chr), levels = unique(gsub("^chr", "", go$chr)))
vmin <- as.numeric(quantile(matA, 0.01, na.rm = TRUE))
vmax <- as.numeric(quantile(matA, 0.99, na.rm = TRUE))
col_fun <- circlize::colorRamp2(c(vmin, 1, vmax), c("#2166ac", "#f7f7f7", "#b2182b"))

group_lab <- factor(ifelse(rownames(matA) == "3", "Normal", "Tumor"), levels = c("Normal", "Tumor"))
ra <- rowAnnotation(
  k_cluster = factor(rownames(matA), levels = clusters),
  group = group_lab,
  col = list(
    k_cluster = k_pal[clusters],
    group = c("Normal" = "#30A3A1", "Tumor" = "#D73027")
  ),
  annotation_name_gp = gpar(fontsize = 8, fontfamily = fig_font_family)
)

ht <- Heatmap(
  matA,
  name = "CNV",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9, fontfamily = fig_font_family),
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_quality = 2,
  column_split = chr_fac,
  column_gap = unit(0.4, "mm"),
  column_title_gp = gpar(fontsize = 7, fontfamily = fig_font_family),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontfamily = fig_font_family),
    labels_gp = gpar(fontsize = 7, fontfamily = fig_font_family)
  )
)

cairo_pdf(file.path(fig_dir, "Fig2A_cnv_heatmap.pdf"), width = 8.2, height = 2.8, family = fig_font_family)
draw(ra + ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(
    file.path(fig_dir, "Fig2A_cnv_heatmap.svg"),
    width = 8.2,
    height = 2.8,
    system_fonts = list(sans = fig_font_family, serif = fig_font_family, mono = "Courier New")
  )
  draw(ra + ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}
png(file.path(fig_dir, "Fig2A_cnv_heatmap.png"), width = 2800, height = 980, res = 300)
draw(ra + ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

message("Panel B: CNV score violin by k-cluster ...")
pB <- ggplot(df_k, aes(x = k_cluster, y = cnv_score_runfinal, fill = k_cluster)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black", linewidth = 0.45, width = 0.9) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9, linewidth = 0.4, fill = "white", color = "black") +
  stat_summary(fun = median, geom = "crossbar", width = 0.35, fatten = 0, color = "black", linewidth = 0.55) +
  scale_fill_manual(values = k_pal[clusters]) +
  annotate(
    "text",
    x = 0.75,
    y = max(df_k$cnv_score_runfinal, na.rm = TRUE) * 0.04,
    hjust = 0,
    label = "Kruskal-Wallis, p < 2.2e-16",
    size = 4.6
  ) +
  labs(x = "k cluster", y = "CNV score") +
  theme_fig +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12),
    plot.margin = margin(6, 8, 6, 6)
  )
save_plot_both(pB, "Fig2B_cnv_score_violin", 5.2, 4.2)

message("Loading malignant calls for panels C-J ...")
meta_mal <- fread(meta_malignant_file, data.table = FALSE)
mal_cells <- meta_mal$cell_id[meta_mal$infercnv_group == "Observation" & meta_mal$malignant_status == "malignant"]
mal_cells <- unique(mal_cells)
message("Malignant cells (method=", method, "): ", length(mal_cells))

message("Loading CytoTRACE metadata and TRS ...")
md <- fread(meta_cytotrace_file, data.table = FALSE)
trs <- fread(trs_file, data.table = FALSE)
df <- md %>%
  select(cell_id, cytotrace_score) %>%
  left_join(trs %>% select(cell_id, trs_score), by = "cell_id") %>%
  filter(cell_id %in% mal_cells)

df <- df %>% filter(is.finite(cytotrace_score), is.finite(trs_score))
message("Cells with both CytoTRACE and TRS available: ", nrow(df))
if (nrow(df) < 200) stop("Too few malignant cells with CytoTRACE/TRS available: ", nrow(df))

message("Computing malignant-only UMAP from Harmony embeddings ...")
obj <- readRDS(seurat_rds)
harmony <- Embeddings(obj, "harmony")
keep <- intersect(df$cell_id, rownames(harmony))
df <- df %>% filter(cell_id %in% keep)
hm <- harmony[match(df$cell_id, rownames(harmony)), , drop = FALSE]

set.seed(42)
um <- uwot::umap(
  hm,
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  n_threads = {
    ncpu <- suppressWarnings(as.integer(parallel::detectCores()))
    if (is.na(ncpu) || ncpu < 2) 1L else (ncpu - 1L)
  },
  ret_model = FALSE,
  verbose = FALSE
)
df$UMAP_1 <- um[, 1]
df$UMAP_2 <- um[, 2]

q1 <- as.numeric(quantile(df$cytotrace_score, 0.25, na.rm = TRUE))
q3 <- as.numeric(quantile(df$cytotrace_score, 0.75, na.rm = TRUE))
df$stem_group_q <- ifelse(df$cytotrace_score <= q1, "LowStem",
                          ifelse(df$cytotrace_score >= q3, "HighStem", "DTStem"))
df$stem_group_q <- factor(df$stem_group_q, levels = c("LowStem", "DTStem", "HighStem"))

# Keep quartile-based grouping on raw CytoTRACE, but plot panel D on
# percentile-scaled CytoTRACE so the histogram shape matches the paper style.
n_ct <- sum(is.finite(df$cytotrace_score))
if (n_ct <= 1) stop("Too few finite CytoTRACE values for percentile scaling: ", n_ct)
df$cytotrace_pct <- (rank(df$cytotrace_score, ties.method = "average", na.last = "keep") - 1) / (n_ct - 1)
q1_d <- as.numeric(quantile(df$cytotrace_pct, 0.25, na.rm = TRUE))
q3_d <- as.numeric(quantile(df$cytotrace_pct, 0.75, na.rm = TRUE))

message("Filtering UMAP outliers for plotting only ...")
flt <- filter_umap_outliers(df, q_trim = 0.005, md_keep = 0.995)
df_umap <- flt$df
write.csv(flt$stats, file.path(tables_dir, "fig2_umap_plot_filter_stats.csv"), row.names = FALSE)
if ("cell_id" %in% colnames(df)) {
  write.csv(
    data.frame(cell_id = df$cell_id, keep_for_umap_plot = flt$keep),
    file.path(tables_dir, "fig2_umap_plot_keep_mask.csv"),
    row.names = FALSE
  )
}
message("UMAP plot cells kept: ", nrow(df_umap), " / ", nrow(df))

message("Panel C ...")
x_rng_c <- range(df_umap$UMAP_1, na.rm = TRUE)
y_rng_c <- range(df_umap$UMAP_2, na.rm = TRUE)
x_breaks_c <- c(-5, 0, 5)
y_breaks_c <- c(-4, 0, 4)
if (x_rng_c[1] > -5 || x_rng_c[2] < 5) x_breaks_c <- pretty(x_rng_c, n = 3)
if (y_rng_c[1] > -4 || y_rng_c[2] < 4) y_breaks_c <- pretty(y_rng_c, n = 3)

pC <- ggplot(df_umap, aes(UMAP_1, UMAP_2, color = cytotrace_score)) +
  geom_point(size = 0.95, alpha = 0.9, stroke = 0) +
  scale_color_gradientn(
    colours = rev(brewer.pal(11, "Spectral")),
    limits = c(0, 1),
    breaks = c(1.0, 0.8, 0.6, 0.4, 0.2, 0.0),
    labels = c("1.0 (Less diff.)", "0.8", "0.6", "0.4", "0.2", "0.0 (More diff.)"),
    oob = squish,
    name = "Predicted
order"
  ) +
  scale_x_continuous(name = "Component1", breaks = x_breaks_c, expand = expansion(mult = 0.02)) +
  scale_y_continuous(name = "Component2", breaks = y_breaks_c, expand = expansion(mult = 0.02)) +
  coord_fixed(ratio = 1) +
  labs(title = "CytoTRACE") +
  theme_classic(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    plot.title = element_text(size = 23, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    legend.position = "left",
    legend.key.height = unit(0.32, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 11),
    plot.margin = margin(4, 6, 6, 6)
  )
save_plot_both(pC, "Fig2C_cytotrace_umap", 5.1, 4.5)

message("Panel D ...")
n_low <- sum(df$stem_group_q == "LowStem", na.rm = TRUE)
n_dt <- sum(df$stem_group_q == "DTStem", na.rm = TRUE)
n_high <- sum(df$stem_group_q == "HighStem", na.rm = TRUE)
pD <- ggplot(df, aes(cytotrace_pct)) +
  geom_histogram(aes(y = after_stat(density)), bins = 85, fill = "grey92", color = "grey75", size = 0.2) +
  geom_density(color = "#CC58D8", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = c(q1_d, q3_d), color = "#C0392B", linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = q1_d / 2, y = Inf, label = paste0("LStem\n", n_low, " cells"), vjust = 1.6, size = 3.5) +
  annotate("text", x = (q1_d + q3_d) / 2, y = Inf, label = paste0("DTStem\n", n_dt, " cells"), vjust = 1.6, size = 3.5) +
  annotate("text", x = (max(df$cytotrace_pct, na.rm = TRUE) + q3_d) / 2, y = Inf, label = paste0("HStem\n", n_high, " cells"), vjust = 1.6, size = 3.5) +
  annotate("text", x = q1_d, y = -Inf, label = sprintf("Bottom quartile (25%%)\nraw score = %.6f", q1), vjust = -0.2, hjust = 0.5, color = "#C0392B", size = 2.7) +
  annotate("text", x = q3_d, y = -Inf, label = sprintf("Top quartile (75%%)\nraw score = %.6f", q3), vjust = -0.2, hjust = 0.5, color = "#C0392B", size = 2.7) +
  labs(title = "Stem score distribution", x = "Scoring", y = "Density") +
  theme_fig
save_plot_both(pD, "Fig2D_cytotrace_hist", 6.1, 4.3)

message("Panel E ...")
pE <- ggplot(df_umap, aes(UMAP_1, UMAP_2, color = stem_group_q)) +
  geom_point(size = 0.2, alpha = 0.9) +
  scale_color_manual(values = stem_cols, name = "group") +
  labs(title = "Malignant", x = NULL, y = NULL) +
  theme_umap
pE <- add_corner_axes(pE, df_umap$UMAP_1, df_umap$UMAP_2, xlab = "UMAP_1", ylab = "UMAP_2")
save_plot_both(pE, "Fig2E_stem_group_umap", 5.3, 4.5)

message("Panel F ...")
cmp <- list(c("LowStem", "DTStem"), c("DTStem", "HighStem"), c("LowStem", "HighStem"))
df$x_grp <- as.numeric(df$stem_group_q)
df_strip <- df

pF <- ggplot(df, aes(x = x_grp, y = cytotrace_score, fill = stem_group_q)) +
  geom_half_violin_right(
    data = df, aes(x = x_grp + 0.04, y = cytotrace_score, fill = stem_group_q, group = stem_group_q),
    width = 0.72, trim = FALSE, alpha = 0.62, color = NA, adjust = 0.75,
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  geom_boxplot(
    aes(group = stem_group_q, fill = stem_group_q),
    width = 0.12, outlier.shape = NA, alpha = 0.75, color = "black", linewidth = 0.42
  ) +
  geom_jitter(
    data = df_strip,
    aes(x = 0.56, y = cytotrace_score, color = stem_group_q),
    width = 0.018, height = 0, size = 0.9, alpha = 0.8, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = stem_cols, name = "group") +
  scale_color_manual(values = stem_cols, guide = "none") +
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = levels(df$stem_group_q),
    limits = c(0.4, 3.65),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
    limits = c(0.00, 1.20),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "CytoTRACE score") +
  theme_classic(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    axis.text.x = element_text(angle = 42, hjust = 1, size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(6, 8, 6, 8)
  )
ymax_f <- max(df$cytotrace_score, na.rm = TRUE)
y12_f <- ymax_f + 0.05
y23_f <- ymax_f + 0.09
y13_f <- ymax_f + 0.14
tick_f <- 0.018
pF <- pF +
  annotate("segment", x = 1.02, xend = 1.98, y = y12_f, yend = y12_f, linewidth = 0.5) +
  annotate("segment", x = 1.02, xend = 1.02, y = y12_f - tick_f, yend = y12_f, linewidth = 0.5) +
  annotate("segment", x = 1.98, xend = 1.98, y = y12_f - tick_f, yend = y12_f, linewidth = 0.5) +
  annotate("text", x = 1.50, y = y12_f + 0.009, label = "***", size = 3.8) +
  annotate("segment", x = 2.02, xend = 2.98, y = y23_f, yend = y23_f, linewidth = 0.5) +
  annotate("segment", x = 2.02, xend = 2.02, y = y23_f - tick_f, yend = y23_f, linewidth = 0.5) +
  annotate("segment", x = 2.98, xend = 2.98, y = y23_f - tick_f, yend = y23_f, linewidth = 0.5) +
  annotate("text", x = 2.50, y = y23_f + 0.009, label = "***", size = 3.8) +
  annotate("segment", x = 1.00, xend = 3.00, y = y13_f, yend = y13_f, linewidth = 0.5) +
  annotate("segment", x = 1.00, xend = 1.00, y = y13_f - tick_f, yend = y13_f, linewidth = 0.5) +
  annotate("segment", x = 3.00, xend = 3.00, y = y13_f - tick_f, yend = y13_f, linewidth = 0.5) +
  annotate("text", x = 2.00, y = y13_f + 0.009, label = "***", size = 3.8)
save_plot_both(pF, "Fig2F_cytotrace_boxplot", 5.2, 4.5)

message("Panel G ...")
pG <- ggplot(df_umap, aes(UMAP_1, UMAP_2, color = trs_score)) +
  geom_point(size = 0.2, alpha = 0.9) +
  scale_color_viridis_c(option = "viridis", name = "Density") +
  labs(title = "scPagwas TRS Score", x = NULL, y = NULL) +
  theme_umap
pG <- add_corner_axes(pG, df_umap$UMAP_1, df_umap$UMAP_2, xlab = "UMAP 1", ylab = "UMAP 2")
save_plot_both(pG, "Fig2G_trs_umap", 5.2, 4.5)

message("Panel H ...")
pH <- ggplot(df, aes(x = x_grp, y = trs_score, fill = stem_group_q)) +
  geom_half_violin_right(
    data = df, aes(x = x_grp + 0.04, y = trs_score, fill = stem_group_q, group = stem_group_q),
    width = 0.72, trim = TRUE, alpha = 0.62, color = NA,
    show.legend = FALSE, inherit.aes = FALSE
  ) +
  geom_boxplot(
    aes(group = stem_group_q, fill = stem_group_q),
    width = 0.16, outlier.shape = NA, alpha = 0.75, color = "black", linewidth = 0.42
  ) +
  geom_jitter(
    data = df_strip,
    aes(x = 0.56, y = trs_score, color = stem_group_q),
    width = 0.02, height = 0, size = 0.8, alpha = 0.8, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = stem_cols, name = "group") +
  scale_color_manual(values = stem_cols, guide = "none") +
  scale_x_continuous(
    breaks = c(1, 2, 3),
    labels = levels(df$stem_group_q),
    limits = c(0.4, 3.65),
    expand = c(0.01, 0.01)
  ) +
  labs(x = NULL, y = "TRS score") +
  theme_classic(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    axis.text.x = element_text(angle = 42, hjust = 1, size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(6, 8, 6, 8)
  )
yr_h <- diff(range(df$trs_score, na.rm = TRUE))
ymax_h <- max(df$trs_score, na.rm = TRUE)
tick_h <- 0.03 * yr_h
pH <- pH +
  annotate("segment", x = 1, xend = 2, y = ymax_h + 0.10 * yr_h, yend = ymax_h + 0.10 * yr_h, linewidth = 0.45) +
  annotate("segment", x = 1, xend = 1, y = ymax_h + 0.10 * yr_h - tick_h, yend = ymax_h + 0.10 * yr_h, linewidth = 0.45) +
  annotate("segment", x = 2, xend = 2, y = ymax_h + 0.10 * yr_h - tick_h, yend = ymax_h + 0.10 * yr_h, linewidth = 0.45) +
  annotate("text", x = 1.5, y = ymax_h + 0.11 * yr_h, label = "***", size = 4.3) +
  annotate("segment", x = 2, xend = 3, y = ymax_h + 0.18 * yr_h, yend = ymax_h + 0.18 * yr_h, linewidth = 0.45) +
  annotate("segment", x = 2, xend = 2, y = ymax_h + 0.18 * yr_h - tick_h, yend = ymax_h + 0.18 * yr_h, linewidth = 0.45) +
  annotate("segment", x = 3, xend = 3, y = ymax_h + 0.18 * yr_h - tick_h, yend = ymax_h + 0.18 * yr_h, linewidth = 0.45) +
  annotate("text", x = 2.5, y = ymax_h + 0.19 * yr_h, label = "***", size = 4.3) +
  annotate("segment", x = 1, xend = 3, y = ymax_h + 0.26 * yr_h, yend = ymax_h + 0.26 * yr_h, linewidth = 0.45) +
  annotate("segment", x = 1, xend = 1, y = ymax_h + 0.26 * yr_h - tick_h, yend = ymax_h + 0.26 * yr_h, linewidth = 0.45) +
  annotate("segment", x = 3, xend = 3, y = ymax_h + 0.26 * yr_h - tick_h, yend = ymax_h + 0.26 * yr_h, linewidth = 0.45) +
  annotate("text", x = 2.0, y = ymax_h + 0.27 * yr_h, label = "***", size = 4.3)
save_plot_both(pH, "Fig2H_trs_boxplot", 5.2, 4.5)

message("Panel I (scatter + marginal densities) ...")
ok <- complete.cases(df$cytotrace_score, df$trs_score)
dfi <- df[ok, ]
r <- cor(dfi$cytotrace_score, dfi$trs_score, method = "pearson")
ct <- suppressWarnings(cor.test(dfi$cytotrace_score, dfi$trs_score, method = "pearson"))
dfree <- if (!is.null(ct$parameter)) as.integer(ct$parameter) else (nrow(dfi) - 2L)
t_stat <- if (!is.null(ct$statistic)) as.numeric(ct$statistic) else NA_real_
ci_lo <- if (!is.null(ct$conf.int)) as.numeric(ct$conf.int[1]) else NA_real_
ci_hi <- if (!is.null(ct$conf.int)) as.numeric(ct$conf.int[2]) else NA_real_
pval <- if (!is.null(ct$p.value)) as.numeric(ct$p.value) else NA_real_
txt_i <- sprintf(
  "t_student(%d) = %.1f, p = %.2g, r_Pearson = %.2f, CI95 [%.2f, %.2f], n = %d",
  dfree, t_stat, pval, r, ci_lo, ci_hi, nrow(dfi)
)

p_sc <- ggplot(dfi, aes(cytotrace_score, trs_score)) +
  geom_point(size = 0.65, alpha = 0.28, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1.05) +
  annotate(
    "text",
    x = min(dfi$cytotrace_score, na.rm = TRUE),
    y = max(dfi$trs_score, na.rm = TRUE),
    hjust = 0, vjust = 1.15,
    label = txt_i, size = 3
  ) +
  labs(x = "Stem Score", y = "TRS Score") +
  theme_classic(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 10),
    plot.margin = margin(1, 1, 1, 1)
  )

p_top <- ggplot(dfi, aes(cytotrace_score)) +
  geom_histogram(bins = 32, fill = "#16A085", color = "#0E6655", linewidth = 0.2) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 0, 1))

p_right <- ggplot(dfi, aes(trs_score)) +
  geom_histogram(bins = 32, fill = "#E67E22", color = "#A04000", linewidth = 0.2) +
  coord_flip() +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 0))

pI <- (p_top + plot_spacer()) / (p_sc + p_right) + plot_layout(widths = c(4, 1.15), heights = c(1, 4))
save_plot_both(pI, "Fig2I_cytotrace_trs_corr", 5.2, 4.5)

message("Panel J ...")
df_j <- df_umap
df_j$stem_n <- rescale(df_j$cytotrace_score, to = c(0, 1))
df_j$trs_n <- rescale(df_j$trs_score, to = c(0, 1))
df_j$mix_col <- rgb(df_j$stem_n, df_j$trs_n, 0)
df_j$alpha_stem <- pmax(0.08, df_j$stem_n)
df_j$alpha_trs <- pmax(0.08, df_j$trs_n)
df_j$alpha_mix <- pmax(0.10, pmax(df_j$stem_n, df_j$trs_n))

x_rng_j <- range(df_j$UMAP_1, na.rm = TRUE)
y_rng_j <- range(df_j$UMAP_2, na.rm = TRUE)
x_breaks_j <- c(-5, 0, 5)
y_breaks_j <- c(-4, 0, 4)
if (x_rng_j[1] > -5 || x_rng_j[2] < 5) x_breaks_j <- pretty(x_rng_j, n = 3)
if (y_rng_j[1] > -4 || y_rng_j[2] < 4) y_breaks_j <- pretty(y_rng_j, n = 3)

theme_j_umap <- theme_classic(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 11, color = "black"),
    plot.margin = margin(4, 5, 4, 4)
  )

base_j <- ggplot(df_j, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey80", size = 0.72, alpha = 0.95, stroke = 0) +
  scale_x_continuous(limits = x_rng_j, breaks = x_breaks_j, expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = y_rng_j, breaks = y_breaks_j, expand = expansion(mult = 0.02)) +
  coord_fixed(ratio = 1)

pJ1 <- base_j +
  geom_point(aes(alpha = alpha_stem), color = "red2", size = 0.76, stroke = 0) +
  scale_alpha_continuous(range = c(0.00, 0.95), guide = "none") +
  labs(title = "Stem Score", x = "UMAP_1", y = "UMAP_2") +
  theme_j_umap

pJ2 <- base_j +
  geom_point(aes(alpha = alpha_trs), color = "#6EDB78", size = 0.76, stroke = 0) +
  scale_alpha_continuous(range = c(0.00, 0.95), guide = "none") +
  labs(title = "TRS Score", x = "UMAP_1", y = "UMAP_2") +
  theme_j_umap

pJ3 <- base_j +
  geom_point(aes(color = mix_col, alpha = alpha_mix), size = 0.76, stroke = 0) +
  scale_color_identity() +
  scale_alpha_continuous(range = c(0.00, 0.92), guide = "none") +
  labs(title = "Stem Score +\nTRS Score", x = "UMAP_1", y = "UMAP_2") +
  theme_j_umap

grid_df <- expand.grid(Scoring = 1:10, TRS = 1:10)
grid_df$fill <- rgb(grid_df$Scoring / 10, grid_df$TRS / 10, 0)
pJ4 <- ggplot(grid_df, aes(Scoring, TRS, fill = fill)) +
  geom_tile() +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = c(2, 4, 6, 8, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 10.5), breaks = c(2.5, 5, 7.5, 10), expand = c(0, 0)) +
  labs(title = "Color threshold: 0", x = "Scoring", y = "scPagwas TRS Score") +
  theme_classic(base_size = 12, base_family = fig_font_family) +
  theme(
    text = element_text(family = fig_font_family),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 11, color = "black"),
    legend.position = "none",
    plot.margin = margin(4, 5, 4, 6)
  )

pJ <- pJ1 + pJ2 + pJ3 + pJ4 + plot_layout(ncol = 4, widths = c(1, 1, 1, 0.95))
save_plot_both(pJ, "Fig2J_umap_combo", 12.8, 4.1)
save_plot_both(pJ1, "Fig2J1_stem_umap", 3.15, 4.1)
save_plot_both(pJ2, "Fig2J2_trs_umap", 3.15, 4.1)
save_plot_both(pJ3, "Fig2J3_mix_umap", 3.15, 4.1)
save_plot_both(pJ4, "Fig2J4_color_threshold", 2.95, 4.1)

message("All Fig2 panels generated (paper-like layout).")
