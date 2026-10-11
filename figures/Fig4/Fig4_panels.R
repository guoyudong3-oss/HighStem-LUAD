#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(WGCNA)
  library(corrplot)
  library(ggplotify)
  library(ggforce)
  library(grid)
})

disableWGCNAThreads()

module_palette <- c(
  turquoise = "#47C9C0",
  green = "#18E318",
  blue = "#1F21FF",
  brown = "#B12A2A",
  yellow = "#F3EA00"
)

target_modules_main <- c("turquoise", "yellow", "green", "brown", "blue")
target_modules_c <- c("turquoise", "green", "blue", "brown", "yellow")

bg_col <- "#FFFFFF"

parse_args <- function(args) {
  run_dir <- NA_character_
  input_rds <- NA_character_
  panel_ref_dir <- NA_character_
  out_prefix <- "Figure4_originalstyle_like"

  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--run-dir" && i < length(args)) {
      run_dir <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--run-dir=")) {
      run_dir <- sub("^--run-dir=", "", a)
      i <- i + 1
      next
    }
    if (a == "--input-rds" && i < length(args)) {
      input_rds <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--input-rds=")) {
      input_rds <- sub("^--input-rds=", "", a)
      i <- i + 1
      next
    }
    if (a == "--panel-ref-dir" && i < length(args)) {
      panel_ref_dir <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--panel-ref-dir=")) {
      panel_ref_dir <- sub("^--panel-ref-dir=", "", a)
      i <- i + 1
      next
    }
    if (a == "--out-prefix" && i < length(args)) {
      out_prefix <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--out-prefix=")) {
      out_prefix <- sub("^--out-prefix=", "", a)
      i <- i + 1
      next
    }
    i <- i + 1
  }

  if (is.na(run_dir) || nchar(run_dir) == 0) {
    stop("Usage: Rscript 02c_make_fig4_original_style_from_run.R --run-dir <path> [--input-rds <seurat_rds>] [--panel-ref-dir <path>] [--out-prefix <name>]")
  }

  list(
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    input_rds = if (!is.na(input_rds)) normalizePath(input_rds, winslash = "/", mustWork = FALSE) else NA_character_,
    panel_ref_dir = if (!is.na(panel_ref_dir)) normalizePath(panel_ref_dir, winslash = "/", mustWork = FALSE) else NA_character_,
    out_prefix = out_prefix
  )
}

guess_input_rds_from_logs <- function(run_dir) {
  logs_dir <- file.path(run_dir, "logs")
  if (!dir.exists(logs_dir)) {
    return(NA_character_)
  }
  log_files <- list.files(logs_dir, pattern = "02b_run_hdwgcna_from_newmal_rds.*\\.log$", full.names = TRUE)
  if (length(log_files) == 0) {
    return(NA_character_)
  }
  log_files <- log_files[order(file.info(log_files)$mtime, decreasing = TRUE)]
  for (lf in log_files) {
    lines <- readLines(lf, warn = FALSE)
    hit <- grep("^Using input rds:", lines, value = TRUE)
    if (length(hit) > 0) {
      p <- trimws(sub("^Using input rds:", "", hit[[1]]))
      if (file.exists(p)) {
        return(normalizePath(p, winslash = "/", mustWork = TRUE))
      }
    }
  }
  NA_character_
}

load_tom_dist <- function(tom_rda) {
  e <- new.env(parent = emptyenv())
  load(tom_rda, envir = e)
  if (exists("consTomDS", envir = e, inherits = FALSE)) {
    obj <- get("consTomDS", envir = e)
    if (inherits(obj, "dist")) {
      return(obj)
    }
    if (is.environment(obj) && exists("consTomDS", envir = obj, inherits = FALSE)) {
      inner <- get("consTomDS", envir = obj)
      if (inherits(inner, "dist")) {
        return(inner)
      }
    }
  }
  stop("Unable to load TOM dist object from: ", tom_rda)
}

get_layout_ratios <- function(panel_ref_dir) {
  out <- list(widths = c(1, 1), heights = c(1, 1, 1.2))
  if (is.na(panel_ref_dir) || nchar(panel_ref_dir) == 0) {
    return(out)
  }
  meta_csv <- file.path(panel_ref_dir, "metadata.csv")
  if (!file.exists(meta_csv)) {
    return(out)
  }
  md <- fread(meta_csv)
  need <- c("panel_id", "w", "h")
  if (!all(need %in% colnames(md))) {
    return(out)
  }
  get_wh <- function(pid) {
    row <- md[panel_id == pid]
    if (nrow(row) == 0) {
      return(c(NA_real_, NA_real_))
    }
    c(as.numeric(row$w[1]), as.numeric(row$h[1]))
  }
  wa <- get_wh("A")[1]; ha <- get_wh("A")[2]
  wb <- get_wh("B")[1]; hb <- get_wh("B")[2]
  wc <- get_wh("C")[1]; hc <- get_wh("C")[2]
  wd <- get_wh("D")[1]; hd <- get_wh("D")[2]
  we <- get_wh("E")[1]; he <- get_wh("E")[2]
  wf <- get_wh("F")[1]; hf <- get_wh("F")[2]
  if (all(is.finite(c(wa, wb, wc, wd, we, wf)))) {
    out$widths <- c(mean(c(wa, wc, we)), mean(c(wb, wd, wf)))
  }
  if (all(is.finite(c(ha, hb, hc, hd, he, hf)))) {
    out$heights <- c(max(ha, hb), max(hc, hd), max(he, hf))
  }
  out
}

get_panel_geometry <- function(panel_ref_dir) {
  fallback <- data.frame(
    panel_id = c("A", "B", "C", "D", "E", "F"),
    x = c(0, 1000, 0, 1000, 0, 1000),
    y = c(0, 0, 760, 760, 1520, 1520),
    w = c(1000, 1000, 1000, 1000, 1000, 1000),
    h = c(740, 740, 740, 740, 860, 860),
    stringsAsFactors = FALSE
  )
  if (is.na(panel_ref_dir) || nchar(panel_ref_dir) == 0) {
    return(fallback)
  }
  meta_csv <- file.path(panel_ref_dir, "metadata.csv")
  if (!file.exists(meta_csv)) {
    return(fallback)
  }
  md <- fread(meta_csv)
  need <- c("panel_id", "x", "y", "w", "h")
  if (!all(need %in% colnames(md))) {
    return(fallback)
  }
  md <- as.data.frame(md[, ..need])
  md <- md[md$panel_id %in% fallback$panel_id, , drop = FALSE]
  if (nrow(md) != 6) {
    return(fallback)
  }
  md$panel_id <- as.character(md$panel_id)
  md$x <- as.numeric(md$x)
  md$y <- as.numeric(md$y)
  md$w <- as.numeric(md$w)
  md$h <- as.numeric(md$h)
  md <- md[match(fallback$panel_id, md$panel_id), , drop = FALSE]
  md
}

scale_to_unit <- function(x, clip = 2) {
  z <- as.numeric(scale(x))
  z[is.na(z)] <- 0
  z <- pmax(pmin(z, clip), -clip) / clip
  z
}

cor_pmat <- function(mat) {
  n <- ncol(mat)
  p_mat <- matrix(NA_real_, n, n)
  diag(p_mat) <- 0
  colnames(p_mat) <- colnames(mat)
  rownames(p_mat) <- colnames(mat)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i < j) {
        ct <- suppressWarnings(cor.test(mat[, i], mat[, j], use = "pairwise.complete.obs"))
        p_mat[i, j] <- ct$p.value
        p_mat[j, i] <- ct$p.value
      }
    }
  }
  p_mat
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
run_dir <- cli$run_dir
if (!dir.exists(run_dir)) {
  stop("Run dir not found: ", run_dir)
}
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
figures_dir <- file.path(run_dir, "figures")
tables_dir <- file.path(run_dir, "tables")
tom_rda <- file.path(run_dir, "tom", "TOM", "HighStem_TOM_TOM.rda")
if (!file.exists(tom_rda)) {
  stop("TOM file not found: ", tom_rda)
}
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- cli$input_rds
if (is.na(input_rds) || !file.exists(input_rds)) {
  guessed <- guess_input_rds_from_logs(run_dir)
  if (is.na(guessed) || !file.exists(guessed)) {
    stop("input_rds not provided/found. Pass --input-rds explicitly.")
  }
  input_rds <- guessed
}
message("Using run dir: ", run_dir)
message("Using input rds: ", input_rds)

power_path <- file.path(tables_dir, "wgcna_power_table.csv")
soft_power_path <- file.path(tables_dir, "wgcna_selected_soft_power.txt")
modules_path <- file.path(tables_dir, "wgcna_modules.csv")
mes_path <- file.path(tables_dir, "wgcna_module_eigengenes_cells.csv")
corr_path <- file.path(tables_dir, "wgcna_module_correlation.csv")
for (p in c(power_path, soft_power_path, modules_path, mes_path)) {
  if (!file.exists(p)) {
    stop("Missing required file: ", p)
  }
}

power_tbl <- fread(power_path)
soft_power <- as.numeric(readLines(soft_power_path, warn = FALSE)[1])
modules <- fread(modules_path)
mes <- fread(mes_path)
corr_tbl <- if (file.exists(corr_path)) fread(corr_path) else NULL
if (!"cell_id" %in% colnames(mes)) {
  setnames(mes, 1, "cell_id")
}
if (!"cell_id" %in% colnames(mes)) {
  stop("wgcna_module_eigengenes_cells.csv missing cell_id column.")
}

for (m in target_modules_main) {
  if (!m %in% colnames(mes)) {
    stop("Missing module eigengene column in MEs table: ", m)
  }
}

seu <- readRDS(input_rds)
if (!inherits(seu, "Seurat")) {
  stop("input_rds is not a Seurat object: ", input_rds)
}
if (!"umap" %in% names(seu@reductions)) {
  stop("Seurat object missing UMAP reduction.")
}
if (!"cytotrace_score" %in% colnames(seu@meta.data)) {
  stop("Seurat metadata missing cytotrace_score.")
}

layout_ratios <- get_layout_ratios(cli$panel_ref_dir)
layout_ratios$widths <- layout_ratios$widths / mean(layout_ratios$widths)
layout_ratios$heights <- layout_ratios$heights / mean(layout_ratios$heights)

theme_base <- theme_gray(base_size = 9) +
  theme(
    plot.background = element_rect(fill = bg_col, color = bg_col),
    panel.background = element_rect(fill = bg_col, color = bg_col),
    legend.background = element_rect(fill = bg_col, color = bg_col),
    legend.key = element_rect(fill = bg_col, color = bg_col),
    strip.background = element_rect(fill = bg_col, color = bg_col)
  )


# ======== Panels drawn at final size (text >= 8 pt) ========
OUT <- "results/reproduce_paper/runs/final_malignant_v3_20261009/W/Fig4/panels/"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
theme_w <- theme_classic(base_size = 8, base_family = "Arial") + theme(text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"))

# ---- A: gene dendrogram and module colours ----
tom_dist <- load_tom_dist(tom_rda); n_tom <- attr(tom_dist, "Size")
if (nrow(modules) > n_tom) modules <- modules[seq_len(n_tom)]
mod_colors <- ifelse(modules$module %in% target_modules_main, modules$module, "grey")
tom_dist_plot <- tom_dist; tom_dist_plot[] <- 1 - tom_dist_plot[]
gene_tree <- fastcluster::hclust(tom_dist_plot, method = "average")
cairo_pdf(paste0(OUT, "Fig4A_dendrogram.pdf"), width = 180 / 25.4, height = 58 / 25.4, pointsize = 8, family = "Arial")
par(mar = c(0.8, 4.0, 1.6, 0.6), cex.axis = 1, cex.lab = 1, cex.main = 1)
WGCNA::plotDendroAndColors(dendro = gene_tree, colors = as.data.frame(mod_colors), groupLabels = "Module", main = "hdWGCNA dendrogram", dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05, guideAll = FALSE, cex.colorLabels = 1, cex.main = 1, cex.axis = 1, cex.lab = 1, lwd = 0.4)
dev.off()

# ---- B: kME curves with top 10 genes per module ----
make_c_panel <- function(mod_name) {
  kcol <- paste0("kME_", mod_name)
  dfm <- as.data.frame(modules) %>% filter(module == mod_name) %>% select(gene_name, all_of(kcol)) %>% arrange(.data[[kcol]])
  plot_df <- data.frame(rank = seq_len(nrow(dfm)), kME = dfm[[kcol]])
  tops <- dfm %>% arrange(desc(.data[[kcol]])) %>% slice_head(n = 10) %>% pull(gene_name)
  y_max <- max(plot_df$kME, na.rm = TRUE)
  text_y <- seq(y_max * 0.80, y_max * 0.10, length.out = length(tops))
  ggplot(plot_df, aes(rank, kME)) + geom_area(fill = module_palette[[mod_name]], color = module_palette[[mod_name]], linewidth = 0.2) +
    annotate("text", x = max(plot_df$rank) * 1.04, y = text_y, label = tops, hjust = 0, size = 2.83, fontface = "italic", color = "black") +
    labs(title = mod_name, x = NULL, y = "kME") + scale_y_continuous(expand = expansion(mult = c(0.01, 0.03))) +
    coord_cartesian(xlim = c(min(plot_df$rank), max(plot_df$rank) * 2.0), clip = "off") + theme_w +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          axis.title.y = element_text(size = 8), axis.text.y = element_text(size = 8), plot.margin = margin(1, 2, 1, 2))
}
pC <- wrap_plots(lapply(target_modules_c, make_c_panel), nrow = 1)
ggsave(paste0(OUT, "Fig4B_kME.pdf"), pC, width = 180 / 25.4, height = 60 / 25.4, device = cairo_pdf)

# ---- C: module activity on UMAP ----
umap <- Embeddings(seu, reduction = "umap")
umap_df <- data.frame(cell_id = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2])
plot_df_d <- merge(umap_df, as.data.frame(mes), by = "cell_id", all = FALSE)
umap_xlim <- as.numeric(quantile(plot_df_d$UMAP_1, probs = c(0.01, 0.99))); umap_ylim <- as.numeric(quantile(plot_df_d$UMAP_2, probs = c(0.02, 0.99)))
x_pad <- diff(umap_xlim) * 0.02; y_pad <- diff(umap_ylim) * 0.02; umap_xlim <- c(umap_xlim[1] - x_pad, umap_xlim[2] + x_pad); umap_ylim <- c(umap_ylim[1] - y_pad, umap_ylim[2] + y_pad)
make_d_panel <- function(mod_name) {
  df <- data.frame(UMAP_1 = plot_df_d$UMAP_1, UMAP_2 = plot_df_d$UMAP_2, score = scale_to_unit(plot_df_d[[mod_name]], clip = 2))
  ggplot(df, aes(UMAP_1, UMAP_2, color = score)) + geom_point(size = 0.3, alpha = 0.82, stroke = 0) +
    scale_color_gradientn(colors = c("#CDCDCD", module_palette[[mod_name]]), limits = c(-1, 1), breaks = c(-1, 1), labels = c("-", "+"), guide = guide_colorbar(title = NULL, barheight = unit(10, "mm"), barwidth = unit(2, "mm"), ticks = FALSE)) +
    labs(title = mod_name) + coord_fixed(xlim = umap_xlim, ylim = umap_ylim, expand = FALSE, clip = "off") + theme_void(base_size = 8) +
    theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold", color = "black", margin = margin(b = 9)), legend.position = "right", legend.text = element_text(size = 8, color = "black"), plot.margin = margin(1, 1, 1, 1))
}
pD <- (make_d_panel("turquoise") + make_d_panel("yellow") + make_d_panel("green")) / (make_d_panel("brown") + make_d_panel("blue") + plot_spacer())
ggsave(paste0(OUT, "Fig4C_module_umap.pdf"), pD, width = 120 / 25.4, height = 86 / 25.4, device = cairo_pdf)

# ---- D: module activity dot plot ----
meta_df <- seu@meta.data; meta_df$cell_id <- rownames(meta_df)
q1 <- as.numeric(quantile(meta_df$cytotrace_score, 0.25, na.rm = TRUE)); q3 <- as.numeric(quantile(meta_df$cytotrace_score, 0.75, na.rm = TRUE))
meta_df$stem_group <- ifelse(meta_df$cytotrace_score >= q3, "HighStem", ifelse(meta_df$cytotrace_score <= q1, "LowStem", "DTStem"))
f_df <- merge(meta_df[, c("cell_id", "stem_group")], as.data.frame(mes)[, c("cell_id", target_modules_main), drop = FALSE], by = "cell_id")
for (m in target_modules_main) f_df[[m]] <- scale_to_unit(f_df[[m]], clip = 2)
f_long <- melt(as.data.table(f_df), id.vars = c("cell_id", "stem_group"), measure.vars = target_modules_main, variable.name = "feature", value.name = "score")
f_sum <- f_long[, .(avg_expr = mean(score, na.rm = TRUE), pct_expr = mean(score > 0, na.rm = TRUE) * 100), by = .(feature, stem_group)]
pct_upper <- max(35, min(75, as.numeric(quantile(f_sum$pct_expr, 0.90, na.rm = TRUE))))
f_sum$pct_expr_plot <- pmin(f_sum$pct_expr, pct_upper)
size_breaks <- c(20, 40, floor(pct_upper)); size_labels <- c("20", "40", paste0(">=", floor(pct_upper)))
f_sum$stem_group <- factor(f_sum$stem_group, levels = c("HighStem", "DTStem", "LowStem")); f_sum$feature <- factor(f_sum$feature, levels = target_modules_main)
hit_rows <- as.character(head(as.data.table(corr_tbl)[module %in% target_modules_main & correlation > 0][order(-correlation)]$module, 2))
y_map <- setNames(seq_along(target_modules_main), target_modules_main)
box_df <- data.frame(xmin = 0.55, xmax = 1.45, ymin = min(y_map[hit_rows]) - 0.45, ymax = max(y_map[hit_rows]) + 0.45)
pF <- ggplot(f_sum, aes(x = stem_group, y = feature, size = pct_expr_plot, color = avg_expr)) + geom_point() +
  geom_rect(data = box_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = NA, color = "red", linetype = "22", linewidth = 0.6) +
  scale_size_continuous(name = "Percent\nexpressed", breaks = size_breaks, labels = size_labels, limits = c(0, pct_upper), range = c(1.3, 5.2)) +
  scale_color_gradient2(name = "Average\nexpression", low = "#169494", mid = "#F1F1F1", high = "#D21F29", midpoint = 0, limits = c(-0.45, 0.45), breaks = c(-0.45, 0, 0.45), oob = scales::squish,
    guide = guide_colorbar(barheight = unit(22, "mm"), barwidth = unit(2.5, "mm"), ticks = TRUE)) +
  labs(x = NULL, y = NULL) + theme_w +
  theme(axis.line = element_line(color = "black", linewidth = 0.4), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 8, color = "black"), legend.position = "right", legend.box = "vertical", legend.spacing.y = unit(1, "mm")) +
  guides(size = guide_legend(order = 1, override.aes = list(color = "black")), color = guide_colorbar(order = 2, barheight = unit(22, "mm"), barwidth = unit(2.5, "mm")))
ggsave(paste0(OUT, "Fig4D_dotplot.pdf"), pF, width = 58 / 25.4, height = 70 / 25.4, device = cairo_pdf)

# ---- soft-threshold and module-correlation panels are drawn by the FigS4 script ----
message("W Fig4 panels done")
