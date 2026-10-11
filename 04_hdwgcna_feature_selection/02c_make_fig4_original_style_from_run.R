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

message("Building panel A ...")
power_metric_plot <- function(metric_col, ylab, hline_fixed = NA_real_, y_limits = NULL, y_breaks = NULL, y_labels = waiver()) {
  df <- power_tbl %>% mutate(value = .data[[metric_col]])
  if (!all(c("Power", "value") %in% colnames(df))) {
    stop("Power table missing required columns for panel A.")
  }
  sel <- df %>% filter(Power == soft_power)
  if (nrow(sel) == 0) {
    stop("Selected soft power not found in power table: ", soft_power)
  }
  hline <- if (is.finite(hline_fixed)) hline_fixed else as.numeric(sel$value[1])
  y_rng <- diff(range(df$value, na.rm = TRUE))
  if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 1
  df <- df %>%
    mutate(
      label_x = case_when(
        Power == soft_power - 3 ~ Power - 0.30,
        Power == soft_power - 2 ~ Power - 0.20,
        Power == soft_power - 1 ~ Power + 0.20,
        TRUE ~ as.numeric(Power)
      ),
      label_y = case_when(
        Power == soft_power ~ value,
        metric_col == "SFT.R.sq" & Power >= 7 ~ pmax(value, 0.88),
        metric_col != "SFT.R.sq" & Power >= 7 ~ pmax(value, 30),
        TRUE ~ value + ifelse((Power %% 2) == 0, y_rng * 0.02, -y_rng * 0.02)
      )
    )

  p <- ggplot(df, aes(x = Power, y = value)) +
    geom_vline(xintercept = soft_power, linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = hline, linetype = "dashed", linewidth = 0.6) +
    geom_text(data = df %>% filter(Power != soft_power), aes(x = label_x, y = label_y, label = Power), size = 1.35, color = "black", check_overlap = TRUE) +
    annotate("point", x = soft_power, y = sel$value[1], size = 3.1, shape = 16, color = "black") +
    annotate("text", x = soft_power, y = sel$value[1], label = as.character(soft_power), color = "white", size = 1.75, fontface = "bold") +
    labs(x = "Soft Power Threshold", y = ylab) +
    scale_x_continuous(limits = c(-0.2, 31), breaks = c(0, 10, 20, 30)) +
    theme_base +
    theme(
      axis.title.x = element_text(size = 5.2, color = "black", margin = margin(t = 2)),
      axis.title.y = element_text(size = 5.2, color = "black", margin = margin(r = 3)),
      axis.text = element_text(size = 4.0, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      panel.grid = element_line(color = "#F2F2F2", linewidth = 0.18),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
      plot.margin = margin(1.5, 2, 1.5, 2)
    )

  if (!is.null(y_limits) || !is.null(y_breaks)) {
    p <- p + scale_y_continuous(limits = y_limits, breaks = y_breaks, labels = y_labels)
  }
  p
}

pA1 <- power_metric_plot(
  "SFT.R.sq", "Scale-free Topology Model",
  hline_fixed = 0.85,
  y_limits = c(-0.05, 1.05),
  y_breaks = seq(0, 1.0, by = 0.2),
  y_labels = function(x) sprintf("%.1f", x)
)
pA2 <- power_metric_plot(
  "mean.k.", "Mean Connectivity",
  y_limits = c(-200, 5300),
  y_breaks = seq(0, 5000, by = 1000),
  y_labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
)
pA3 <- power_metric_plot(
  "median.k.", "Median Connectivity",
  y_limits = c(-200, 5300),
  y_breaks = seq(0, 5000, by = 1000),
  y_labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
)
pA4 <- power_metric_plot(
  "max.k.", "Max Connectivity",
  y_limits = c(-200, 6000),
  y_breaks = seq(0, 6000, by = 2000),
  y_labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
)
pA1 <- pA1 + theme(plot.margin = margin(1.5, 2.0, 1.2, 13))
pA2 <- pA2 + theme(plot.margin = margin(1.5, 2.0, 1.2, 6))
pA3 <- pA3 + theme(plot.margin = margin(1.2, 2.0, 1.5, 13))
pA4 <- pA4 + theme(plot.margin = margin(1.2, 2.0, 1.5, 6))
pA <- (pA1 + pA2) / (pA3 + pA4)

message("Building panel B ...")
tom_dist <- load_tom_dist(tom_rda)
n_tom <- attr(tom_dist, "Size")
if (!is.numeric(n_tom) || is.na(n_tom)) {
  stop("Invalid TOM dist size.")
}
if (nrow(modules) < n_tom) {
  stop("wgcna_modules.csv has fewer rows than TOM size.")
}
if (nrow(modules) > n_tom) {
  modules <- modules[seq_len(n_tom)]
}

mod_colors <- ifelse(modules$module %in% target_modules_main, modules$module, "grey")
# Match original Figure4B axis scale: plot on (1 - TOM distance) to place heights in ~0.65-1 range.
tom_dist_plot <- tom_dist
tom_dist_plot[] <- 1 - tom_dist_plot[]
gene_tree <- fastcluster::hclust(tom_dist_plot, method = "average")
draw_panel_B_base <- function() {
  par(bg = bg_col, mar = c(4.2, 4.4, 2.5, 1.2), las = 0, cex.axis = 1.0, cex.lab = 1.6, cex.main = 1.6)
  WGCNA::plotDendroAndColors(
    dendro = gene_tree,
    colors = as.data.frame(mod_colors),
    groupLabels = "Module colors",
    main = "hdWGCNA Dendrogram",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    guideAll = FALSE,
    cex.rowText = 1.5,
    cex.colorLabels = 1.2,
    lwd = 0.6
  )
}
save_panel_B_base <- function(out_file, width_px, height_px, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    # Keep classic base-graphics sizing semantics explicit.
    png(
      filename = out_file,
      width = width_px,
      height = height_px,
      units = "px",
      res = 72,
      pointsize = 12,
      bg = bg_col
    )
  } else {
    # Match PNG visual scale (72 px ~= 1 inch in default base graphics).
    cairo_pdf(
      filename = out_file,
      width = width_px / 72,
      height = height_px / 72,
      pointsize = 12,
      bg = bg_col
    )
  }
  on.exit(dev.off(), add = TRUE)
  draw_panel_B_base()
}
pB <- ggplot() + theme_void()

message("Building panel C ...")
make_c_panel <- function(mod_name) {
  kcol <- paste0("kME_", mod_name)
  if (!kcol %in% colnames(modules)) {
    return(ggplot() + theme_void() + ggtitle(mod_name))
  }
  dfm <- as.data.frame(modules) %>%
    filter(module == mod_name) %>%
    select(gene_name, all_of(kcol))
  if (nrow(dfm) == 0) {
    return(ggplot() + theme_void() + ggtitle(mod_name))
  }
  dfm <- dfm %>% arrange(.data[[kcol]])
  plot_df <- data.frame(rank = seq_len(nrow(dfm)), kME = dfm[[kcol]])
  tops <- dfm %>% arrange(desc(.data[[kcol]])) %>% slice_head(n = min(10, nrow(dfm))) %>% pull(gene_name)
  y_min <- min(plot_df$kME, na.rm = TRUE)
  y_max <- max(plot_df$kME, na.rm = TRUE)
  text_y <- seq(y_max * 0.68, y_max * 0.18, length.out = length(tops))
  if (!all(is.finite(text_y)) || y_max <= 0) {
    text_y <- seq(y_max, y_min, length.out = length(tops))
  }
  x_right <- max(plot_df$rank) * 1.30
  text_x <- max(plot_df$rank) * 1.03
  ggplot(plot_df, aes(rank, kME)) +
    geom_area(fill = module_palette[[mod_name]], color = module_palette[[mod_name]], alpha = 0.98, linewidth = 0.25) +
    annotate(
      "text",
      x = text_x,
      y = text_y,
      label = tops,
      hjust = 0,
      size = 1.55,
      fontface = "italic",
      color = "black",
      lineheight = 0.86
    ) +
    labs(title = mod_name, x = NULL, y = "kME") +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.03))) +
    coord_cartesian(xlim = c(min(plot_df$rank), x_right), clip = "off") +
    theme_base +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 10.5, hjust = 0.5, color = "black"),
      axis.title.y = element_text(size = 9.0, color = "black"),
      axis.text.y = element_text(size = 7.4, color = "black"),
      axis.ticks.y = element_line(color = "black", linewidth = 0.35),
      plot.margin = margin(1, 14, 1, 2)
    )
}

pC_list <- lapply(target_modules_c, make_c_panel)
pC <- wrap_plots(pC_list, nrow = 1)

message("Building panel D ...")
umap <- Embeddings(seu, reduction = "umap")
umap_df <- data.frame(cell_id = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], stringsAsFactors = FALSE)
plot_df_d <- merge(umap_df, as.data.frame(mes), by = "cell_id", all = FALSE)
if (nrow(plot_df_d) == 0) {
  stop("No overlapping cells between UMAP and module eigengene table.")
}
# Use robust limits for display so a few outliers do not distort panel proportions.
umap_xlim <- as.numeric(quantile(plot_df_d$UMAP_1, probs = c(0.01, 0.99), na.rm = TRUE))
umap_ylim <- as.numeric(quantile(plot_df_d$UMAP_2, probs = c(0.02, 0.99), na.rm = TRUE))
x_pad <- diff(umap_xlim) * 0.02
y_pad <- diff(umap_ylim) * 0.02
umap_xlim <- c(umap_xlim[1] - x_pad, umap_xlim[2] + x_pad)
umap_ylim <- c(umap_ylim[1] - y_pad, umap_ylim[2] + y_pad)

make_d_panel <- function(mod_name) {
  val <- scale_to_unit(plot_df_d[[mod_name]], clip = 2)
  df <- data.frame(
    UMAP_1 = plot_df_d$UMAP_1,
    UMAP_2 = plot_df_d$UMAP_2,
    score = val
  )
  ggplot(df, aes(UMAP_1, UMAP_2, color = score)) +
    geom_point(size = 0.34, alpha = 0.82) +
    scale_color_gradientn(
      colors = c("#CDCDCD", module_palette[[mod_name]]),
      limits = c(-1, 1),
      breaks = c(-1, 1),
      labels = c("-", "+"),
      guide = guide_colorbar(
        title = NULL,
        barheight = unit(1.1, "cm"),
        barwidth = unit(0.18, "cm"),
        ticks = FALSE
      )
    ) +
    labs(title = mod_name) +
    coord_fixed(xlim = umap_xlim, ylim = umap_ylim, expand = FALSE, clip = "off") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = bg_col, color = bg_col),
      panel.background = element_rect(fill = bg_col, color = bg_col),
      plot.title = element_text(size = 10, hjust = 0.5, color = "black"),
      legend.position = "right",
      legend.text = element_text(size = 7, color = "black"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

pD_turquoise <- make_d_panel("turquoise")
pD_yellow <- make_d_panel("yellow")
pD_green <- make_d_panel("green")
pD_brown <- make_d_panel("brown")
pD_blue <- make_d_panel("blue")
pD <- (pD_turquoise + pD_yellow + pD_green) / (pD_brown + pD_blue + plot_spacer())

message("Building panel E ...")
me_sub <- as.matrix(as.data.frame(mes)[, target_modules_main, drop = FALSE])
cor_mat <- cor(me_sub, use = "pairwise.complete.obs")
p_mat <- cor_pmat(me_sub)

mods <- target_modules_main
n_mod <- length(mods)
grid_df <- expand.grid(row = seq_len(n_mod), col = seq_len(n_mod))
grid_df <- grid_df[grid_df$col >= grid_df$row, , drop = FALSE]
grid_df$mod_row <- mods[grid_df$row]
grid_df$mod_col <- mods[grid_df$col]
grid_df$r <- mapply(function(i, j) cor_mat[i, j], grid_df$mod_row, grid_df$mod_col)
grid_df$p <- mapply(function(i, j) p_mat[i, j], grid_df$mod_row, grid_df$mod_col)
grid_df$x <- grid_df$col
grid_df$y <- n_mod - grid_df$row + 1
grid_df$is_diag <- grid_df$row == grid_df$col
grid_df$absr <- pmin(abs(grid_df$r), 0.999)
grid_df$angle <- ifelse(grid_df$r >= 0, pi / 4, -pi / 4)
grid_df$b <- 0.42 * sqrt(1 - grid_df$absr)

diag_df <- grid_df[grid_df$is_diag, , drop = FALSE]
offdiag_df <- grid_df[!grid_df$is_diag, , drop = FALSE]
insig_df <- offdiag_df[offdiag_df$p >= 0.05, , drop = FALSE]

row_label_df <- data.frame(
  label = mods,
  x = 0.50,
  y = rev(seq_len(n_mod)),
  stringsAsFactors = FALSE
)
col_label_df <- data.frame(
  label = mods,
  x = seq_len(n_mod),
  y = n_mod + 0.62,
  stringsAsFactors = FALSE
)

pE <- ggplot() +
  geom_rect(
    data = grid_df,
    aes(xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.5, ymax = y + 0.5),
    fill = NA, color = "#B3B3B3", linewidth = 0.85
  ) +
  ggforce::geom_ellipse(
    data = offdiag_df,
    aes(x0 = x, y0 = y, a = 0.43, b = b, angle = angle, fill = r),
    color = NA, alpha = 0.95
  ) +
  geom_segment(
    data = diag_df,
    aes(x = x - 0.40, y = y - 0.40, xend = x + 0.40, yend = y + 0.40),
    color = "#B742FF", linewidth = 1.05
  ) +
  geom_text(
    data = insig_df,
    aes(x = x, y = y),
    label = "×", size = 8.6, color = "black"
  ) +
  geom_text(
    data = row_label_df,
    aes(x = x, y = y, label = label),
    hjust = 1, vjust = 0.5, size = 2.9, color = "black"
  ) +
  geom_text(
    data = col_label_df,
    aes(x = x, y = y, label = label),
    angle = 45, hjust = 0, vjust = 0.5, size = 2.9, color = "black"
  ) +
  scale_fill_gradientn(
    colors = colorRampPalette(c("#2E8B57", "#F0F0F0", "#B742FF"))(200),
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.2),
    name = NULL,
    guide = guide_colorbar(
      barheight = unit(1.7, "in"),
      barwidth = unit(0.08, "in"),
      frame.colour = "black",
      ticks.colour = "#E6E6E6"
    )
  ) +
  coord_equal(xlim = c(0.20, n_mod + 0.76), ylim = c(0.36, n_mod + 0.98), expand = FALSE, clip = "off") +
  theme_void() +
  theme(
    plot.background = element_rect(fill = bg_col, color = bg_col),
    panel.background = element_rect(fill = bg_col, color = bg_col),
    legend.position = "right",
    legend.text = element_text(size = 7, color = "black"),
    plot.margin = margin(18, 6, 4, 26)
  )

message("Building panel F ...")
meta_df <- seu@meta.data
meta_df$cell_id <- rownames(meta_df)
q1 <- as.numeric(quantile(meta_df$cytotrace_score, 0.25, na.rm = TRUE))
q3 <- as.numeric(quantile(meta_df$cytotrace_score, 0.75, na.rm = TRUE))
meta_df$stem_group <- ifelse(meta_df$cytotrace_score >= q3, "HighStem", ifelse(meta_df$cytotrace_score <= q1, "LowStem", "DTStem"))

f_df <- merge(meta_df[, c("cell_id", "stem_group")], as.data.frame(mes)[, c("cell_id", target_modules_main), drop = FALSE], by = "cell_id", all = FALSE)
for (m in target_modules_main) {
  f_df[[m]] <- scale_to_unit(f_df[[m]], clip = 2)
}

f_long <- melt(
  as.data.table(f_df),
  id.vars = c("cell_id", "stem_group"),
  measure.vars = target_modules_main,
  variable.name = "feature",
  value.name = "score"
)
f_sum <- f_long[, .(
  avg_expr = mean(score, na.rm = TRUE),
  pct_expr = mean(score > 0, na.rm = TRUE) * 100
), by = .(feature, stem_group)]
# Adapt bubble sizing to this run's distribution (avoid oversized bubbles).
pct_upper <- as.numeric(quantile(f_sum$pct_expr, 0.90, na.rm = TRUE))
pct_upper <- max(35, min(75, pct_upper))
f_sum$pct_expr_plot <- pmin(f_sum$pct_expr, pct_upper)
size_breaks <- pretty(c(0, pct_upper), n = 4)
size_breaks <- size_breaks[size_breaks > 0]
if (length(size_breaks) > 3) {
  size_breaks <- tail(size_breaks, 3)
}

f_sum$stem_group <- factor(f_sum$stem_group, levels = c("HighStem", "DTStem", "LowStem"))
f_sum$feature <- factor(f_sum$feature, levels = c("turquoise", "yellow", "green", "brown", "blue"))

stem_levels <- c("HighStem", "DTStem", "LowStem")
feature_levels <- c("turquoise", "yellow", "green", "brown", "blue")
hit_rows <- character(0)
if (!is.null(corr_tbl) && all(c("module", "correlation") %in% colnames(corr_tbl))) {
  corr_hits <- as.data.table(corr_tbl)[
    module %in% feature_levels & is.finite(correlation) & correlation > 0
  ][order(-correlation)]
  if (nrow(corr_hits) > 0) {
    hit_rows <- as.character(head(corr_hits$module, 2))
  }
}
if (length(hit_rows) == 0) {
  fs_wide <- dcast(as.data.table(f_sum), feature ~ stem_group, value.var = "avg_expr")
  hit_rows <- fs_wide[HighStem > pmax(DTStem, LowStem) & HighStem > 0, as.character(feature)]
}
if (length(hit_rows) > 0) {
  hit_cols <- "HighStem"
  x_map <- setNames(seq_along(stem_levels), stem_levels)
  y_map <- setNames(seq_along(feature_levels), feature_levels)
  box_df <- data.frame(
    xmin = min(x_map[hit_cols]) - 0.45,
    xmax = max(x_map[hit_cols]) + 0.45,
    ymin = min(y_map[hit_rows]) - 0.45,
    ymax = max(y_map[hit_rows]) + 0.45
  )
} else {
  box_df <- NULL
}

pF <- ggplot(f_sum, aes(x = stem_group, y = feature, size = pct_expr_plot, color = avg_expr)) +
  geom_point() +
  {
    if (!is.null(box_df)) {
      geom_rect(
        data = box_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = "red",
        linetype = "dashed",
        linewidth = 1.1
      )
    } else {
      NULL
    }
  } +
  scale_size_continuous(
    name = "Percent Expressed",
    breaks = size_breaks,
    limits = c(0, pct_upper),
    range = c(1.3, 5.4)
  ) +
  scale_color_gradient2(
    name = "Average Expression",
    low = "#169494",
    mid = "#F1F1F1",
    high = "#D21F29",
    midpoint = 0,
    limits = c(-0.45, 0.45),
    breaks = c(-0.45, -0.225, 0, 0.225, 0.45),
    oob = scales::squish,
    guide = guide_colorbar(
      barheight = unit(1.65, "in"),
      barwidth = unit(0.16, "in"),
      ticks = TRUE
    )
  ) +
  labs(x = NULL, y = "Features") +
  theme_classic(base_size = 10) +
  theme(
    plot.background = element_rect(fill = bg_col, color = bg_col),
    panel.background = element_rect(fill = bg_col, color = bg_col),
    axis.line = element_line(color = "black", linewidth = 0.9),
    axis.title = element_text(size = 11.5, color = "black"),
    axis.text.x = element_text(size = 7.0, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10.5, color = "black"),
    legend.title = element_text(size = 11.5, color = "black"),
    legend.text = element_text(size = 7.2, color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.14, "cm")
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(color = "black", alpha = 1)),
    color = guide_colorbar(order = 2)
  )

message("Composing full Figure4 layout ...")
panel_files <- list(
  A = pA, B = pB, C = pC, D = pD, E = pE, F = pF
)

panel_geo <- get_panel_geometry(cli$panel_ref_dir)
if (!is.na(cli$panel_ref_dir) && nchar(cli$panel_ref_dir) > 0 && file.exists(file.path(cli$panel_ref_dir, "metadata.csv"))) {
  message("Using panel reference geometry from: ", cli$panel_ref_dir)
} else {
  message("Using fallback panel geometry (no valid --panel-ref-dir metadata.csv found).")
}
for (i in seq_len(nrow(panel_geo))) {
  message(sprintf("Panel %s geometry: x=%d y=%d w=%d h=%d",
                  panel_geo$panel_id[i], as.integer(panel_geo$x[i]), as.integer(panel_geo$y[i]),
                  as.integer(panel_geo$w[i]), as.integer(panel_geo$h[i])))
}
canvas_w <- max(panel_geo$x + panel_geo$w)
canvas_h <- max(panel_geo$y + panel_geo$h)

panel_png_paths <- setNames(
  file.path(figures_dir, paste0(cli$out_prefix, "_", names(panel_files), ".png")),
  names(panel_files)
)
panel_pdf_paths <- setNames(
  file.path(figures_dir, paste0(cli$out_prefix, "_", names(panel_files), ".pdf")),
  names(panel_files)
)

for (nm in names(panel_files)) {
  g <- panel_geo[panel_geo$panel_id == nm, , drop = FALSE]
  w_in <- as.numeric(g$w) / 300
  h_in <- as.numeric(g$h) / 300
  if (nm == "B") {
    save_panel_B_base(panel_png_paths[[nm]], as.numeric(g$w), as.numeric(g$h), device = "png")
    save_panel_B_base(panel_pdf_paths[[nm]], as.numeric(g$w), as.numeric(g$h), device = "pdf")
    next
  }
  ggsave(panel_png_paths[[nm]], panel_files[[nm]], width = w_in, height = h_in, units = "in", dpi = 300, bg = bg_col, limitsize = FALSE)
  ggsave(panel_pdf_paths[[nm]], panel_files[[nm]], width = w_in, height = h_in, units = "in", dpi = 300, bg = bg_col, limitsize = FALSE)
}

compose_canvas <- function(out_file, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(filename = out_file, width = canvas_w, height = canvas_h, bg = bg_col)
  } else {
    cairo_pdf(filename = out_file, width = canvas_w / 300, height = canvas_h / 300, bg = bg_col)
  }
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  grid.rect(gp = gpar(fill = bg_col, col = bg_col))
  for (i in seq_len(nrow(panel_geo))) {
    pid <- panel_geo$panel_id[i]
    img <- png::readPNG(panel_png_paths[[pid]])
    x_center <- (panel_geo$x[i] + panel_geo$w[i] / 2) / canvas_w
    y_center <- 1 - (panel_geo$y[i] + panel_geo$h[i] / 2) / canvas_h
    grid.raster(
      img,
      x = unit(x_center, "npc"),
      y = unit(y_center, "npc"),
      width = unit(panel_geo$w[i] / canvas_w, "npc"),
      height = unit(panel_geo$h[i] / canvas_h, "npc"),
      interpolate = FALSE
    )
    grid.text(
      pid,
      x = unit((panel_geo$x[i] + 12) / canvas_w, "npc"),
      y = unit(1 - (panel_geo$y[i] + 8) / canvas_h, "npc"),
      just = c("left", "top"),
      gp = gpar(col = "black", fontsize = 15, fontface = "bold", fontfamily = "serif")
    )
  }
}

out_png <- file.path(figures_dir, paste0(cli$out_prefix, ".png"))
out_pdf <- file.path(figures_dir, paste0(cli$out_prefix, ".pdf"))
compose_canvas(out_png, device = "png")
compose_canvas(out_pdf, device = "pdf")

# Also export an editable vector full figure with exact panel geometry (for Illustrator refinement).
panel_grobs <- list(
  A = patchwork::patchworkGrob(pA),
  # Keep panel B editable in vector output by drawing the base dendrogram as a grob.
  B = ggplotify::as.grob(draw_panel_B_base),
  C = patchwork::patchworkGrob(pC),
  D = patchwork::patchworkGrob(pD),
  E = ggplot2::ggplotGrob(pE),
  F = ggplot2::ggplotGrob(pF)
)

compose_vector_canvas <- function(out_file) {
  cairo_pdf(filename = out_file, width = canvas_w / 300, height = canvas_h / 300, bg = bg_col)
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  grid.rect(gp = gpar(fill = bg_col, col = bg_col))
  for (i in seq_len(nrow(panel_geo))) {
    pid <- panel_geo$panel_id[i]
    x_center <- (panel_geo$x[i] + panel_geo$w[i] / 2) / canvas_w
    y_center <- 1 - (panel_geo$y[i] + panel_geo$h[i] / 2) / canvas_h
    vp <- viewport(
      x = unit(x_center, "npc"),
      y = unit(y_center, "npc"),
      width = unit(panel_geo$w[i] / canvas_w, "npc"),
      height = unit(panel_geo$h[i] / canvas_h, "npc"),
      just = c("center", "center"),
      clip = "off"
    )
    pushViewport(vp)
    grid.draw(panel_grobs[[pid]])
    popViewport()
    grid.text(
      pid,
      x = unit((panel_geo$x[i] + 12) / canvas_w, "npc"),
      y = unit(1 - (panel_geo$y[i] + 8) / canvas_h, "npc"),
      just = c("left", "top"),
      gp = gpar(col = "black", fontsize = 15, fontface = "bold", fontfamily = "serif")
    )
  }
}

out_edit_pdf <- file.path(figures_dir, paste0(cli$out_prefix, "_editable_vector.pdf"))
compose_vector_canvas(out_edit_pdf)
out_edit_ai <- file.path(figures_dir, paste0(cli$out_prefix, "_editable_vector.ai"))
file.copy(out_edit_pdf, out_edit_ai, overwrite = TRUE)

message("Saved full Figure4 (style-matched): ", out_pdf)
message("Saved full Figure4 (style-matched): ", out_png)
message("Saved full Figure4 (editable vector): ", out_edit_pdf)
message("Saved full Figure4 (editable AI copy): ", out_edit_ai)
message("Saved A-F panel PNG/PDF with prefix: ", cli$out_prefix)
message("Done.")
