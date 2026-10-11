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


OUT <- "results/reproduce_paper/runs/final_malignant_v3_20261009/W/FigS4/"
theme_w <- theme_bw(base_size = 8, base_family = "Arial") + theme(text = element_text(size = 8, color = "black"), axis.text = element_text(size = 8, color = "black"), panel.grid = element_line(color = "#F2F2F2", linewidth = 0.18), panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5), plot.margin = margin(2, 4, 2, 2))
power_metric_plot <- function(metric_col, ylab, hline_fixed = NA_real_, y_limits = NULL, y_breaks = NULL, y_labels = waiver()) {
  df <- power_tbl %>% mutate(value = .data[[metric_col]]); sel <- df %>% filter(Power == soft_power)
  hline <- if (is.finite(hline_fixed)) hline_fixed else as.numeric(sel$value[1])
  p <- ggplot(df, aes(x = Power, y = value)) + geom_vline(xintercept = soft_power, linetype = "dashed", linewidth = 0.4) + geom_hline(yintercept = hline, linetype = "dashed", linewidth = 0.4) +
    geom_point(data = df %>% filter(Power != soft_power), size = 0.5, color = "grey35") + geom_text(data = df %>% filter(Power != soft_power, Power %in% (if (metric_col == "SFT.R.sq") c(1:9, 12, 16, 20, 24, 28) else c(1:6, 12, 16, 20, 24, 28))) %>% mutate(lx = Power + ifelse(Power == 9, -1.3, ifelse(Power == 12, 1.5, 0)), ly = value + ifelse(Power == 9, diff(range(value)) * 0.07, diff(range(value)) * 0.045)), aes(x = lx, y = ly, label = Power), size = 2.83, color = "black") +
    annotate("point", x = soft_power, y = sel$value[1], size = 4.6, shape = 16, color = "black") + annotate("text", x = soft_power, y = sel$value[1], label = as.character(soft_power), color = "white", size = 2.83, fontface = "bold") +
    labs(x = "Soft power threshold", y = ylab) + scale_x_continuous(limits = c(-0.2, 31), breaks = c(0, 10, 20, 30)) + theme_w
  if (!is.null(y_limits) || !is.null(y_breaks)) p <- p + scale_y_continuous(limits = y_limits, breaks = y_breaks, labels = y_labels)
  p
}
cm <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
pA <- (power_metric_plot("SFT.R.sq", "Scale-free topology fit (R²)", 0.85, c(-0.05, 1.05), seq(0, 1, 0.2), function(x) sprintf("%.1f", x)) + power_metric_plot("mean.k.", "Mean connectivity", NA_real_, c(-200, 5300), seq(0, 5000, 1000), cm)) /
      (power_metric_plot("median.k.", "Median connectivity", NA_real_, c(-200, 5300), seq(0, 5000, 1000), cm) + power_metric_plot("max.k.", "Max connectivity", NA_real_, c(-200, 6000), seq(0, 6000, 2000), cm))
ggsave(paste0(OUT, "FigS4A_soft_power_v2.pdf"), pA, width = 130 / 25.4, height = 100 / 25.4, device = cairo_pdf)
# module correlation matrix (ellipses)
me_sub <- as.matrix(as.data.frame(mes)[, target_modules_main, drop = FALSE]); cor_mat <- cor(me_sub, use = "pairwise.complete.obs"); p_mat <- cor_pmat(me_sub)
mods <- target_modules_main; n_mod <- length(mods)
grid_df <- expand.grid(row = seq_len(n_mod), col = seq_len(n_mod)); grid_df <- grid_df[grid_df$col >= grid_df$row, , drop = FALSE]
grid_df$mod_row <- mods[grid_df$row]; grid_df$mod_col <- mods[grid_df$col]
grid_df$r <- mapply(function(i, j) cor_mat[i, j], grid_df$mod_row, grid_df$mod_col); grid_df$p <- mapply(function(i, j) p_mat[i, j], grid_df$mod_row, grid_df$mod_col)
grid_df$x <- grid_df$col; grid_df$y <- n_mod - grid_df$row + 1; grid_df$is_diag <- grid_df$row == grid_df$col; grid_df$absr <- pmin(abs(grid_df$r), 0.999)
grid_df$angle <- ifelse(grid_df$r >= 0, pi / 4, -pi / 4); grid_df$b <- 0.42 * sqrt(1 - grid_df$absr)
diag_df <- grid_df[grid_df$is_diag, , drop = FALSE]; offdiag_df <- grid_df[!grid_df$is_diag, , drop = FALSE]; insig_df <- offdiag_df[offdiag_df$p >= 0.05, , drop = FALSE]
row_label_df <- data.frame(label = mods, x = 0.50, y = rev(seq_len(n_mod))); col_label_df <- data.frame(label = mods, x = seq_len(n_mod), y = n_mod + 0.62)
pE <- ggplot() + geom_rect(data = grid_df, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5), fill = NA, color = "#B3B3B3", linewidth = 0.5) +
  ggforce::geom_ellipse(data = offdiag_df, aes(x0 = x, y0 = y, a = 0.43, b = b, angle = angle, fill = r), color = NA) +
  geom_segment(data = diag_df, aes(x = x - 0.40, y = y - 0.40, xend = x + 0.40, yend = y + 0.40), color = "#B742FF", linewidth = 0.7) +
  geom_text(data = insig_df, aes(x = x, y = y), label = "x", size = 4, color = "black") +
  geom_text(data = row_label_df, aes(x = x, y = y, label = label), hjust = 1, vjust = 0.5, size = 2.83, color = "black") +
  geom_text(data = col_label_df, aes(x = x, y = y, label = label), angle = 45, hjust = 0, vjust = 0.5, size = 2.83, color = "black") +
  scale_fill_gradientn(colors = colorRampPalette(c("#2E8B57", "#F0F0F0", "#B742FF"))(200), limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), name = "Correlation", guide = guide_colorbar(barheight = unit(30, "mm"), barwidth = unit(2.5, "mm"), frame.colour = "black", ticks.colour = "#E6E6E6")) +
  coord_equal(xlim = c(-0.7, n_mod + 1.2), ylim = c(0.36, n_mod + 1.9), expand = FALSE, clip = "off") + theme_void(base_size = 8) +
  theme(legend.position = "right", legend.text = element_text(size = 8, color = "black"), legend.title = element_text(size = 8, color = "black"), plot.margin = margin(2, 2, 2, 2))
ggsave(paste0(OUT, "FigS4E_module_corr_v2.pdf"), pE, width = 90 / 25.4, height = 74 / 25.4, device = cairo_pdf)
message("W FigS4 done")
