#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
  library(png)
  library(jsonlite)
  library(scales)
})

get_arg <- function(args, key, default = NA_character_) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", key, "="), "", hit[[1]]))
  idx <- which(args == key)
  if (length(idx) > 0 && idx[[1]] < length(args)) return(args[idx[[1]] + 1])
  default
}

find_sample_dir <- function(root_dir, sample_name) {
  direct <- file.path(root_dir, "data", "spatial", "GSE307534", sample_name)
  if (dir.exists(direct)) return(normalizePath(direct, winslash = "/", mustWork = TRUE))

  cohort_root <- file.path(root_dir, "data", "spatial", "GSE307534", "Full_Cohort")
  cand <- Sys.glob(file.path(cohort_root, "*", sample_name))
  cand <- cand[dir.exists(cand)]
  if (length(cand) > 0) {
    return(normalizePath(cand[[1]], winslash = "/", mustWork = TRUE))
  }

  stop("Cannot find sample directory for ", sample_name)
}

read_positions <- function(sample_dir) {
  pos_path <- file.path(sample_dir, "spatial", "tissue_positions.csv")
  if (!file.exists(pos_path)) {
    pos_path <- file.path(sample_dir, "spatial", "tissue_positions_list.csv")
  }
  if (!file.exists(pos_path)) stop("Missing tissue positions file in ", sample_dir)

  dt <- fread(pos_path)
  if (all(c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres") %in% names(dt))) {
    return(dt)
  }

  dt <- fread(pos_path, header = FALSE)
  if (ncol(dt) < 6) stop("Unexpected tissue positions format: ", pos_path)
  setnames(dt, 1:6, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"))
  dt
}

make_rgba <- function(img, alpha = 0.28) {
  d <- dim(img)
  if (length(d) == 2) {
    arr <- array(1, dim = c(d[1], d[2], 4))
    arr[, , 1] <- img
    arr[, , 2] <- img
    arr[, , 3] <- img
    arr[, , 4] <- alpha
    return(arr)
  }
  if (d[3] == 4) {
    out <- img
    out[, , 4] <- alpha
    return(out)
  }
  arr <- array(1, dim = c(d[1], d[2], 4))
  arr[, , 1:3] <- img[, , 1:3]
  arr[, , 4] <- alpha
  arr
}

hex_vertices_dt <- function(df, radius) {
  ang <- seq(0, 2 * pi, length.out = 7)[1:6]
  pieces <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    pieces[[i]] <- data.table(
      barcode = df$barcode[[i]],
      Score = df$Score[[i]],
      x = df$x[[i]] + radius * cos(ang),
      y = df$y[[i]] + radius * sin(ang),
      ord = seq_along(ang)
    )
  }
  rbindlist(pieces)
}

pretty_breaks_02 <- function(v) {
  q <- quantile(v, probs = c(0.02, 0.98), na.rm = TRUE)
  lo <- floor(q[[1]] / 0.2) * 0.2
  hi <- ceiling(q[[2]] / 0.2) * 0.2
  if (!is.finite(lo)) lo <- -0.2
  if (!is.finite(hi)) hi <- 0.4
  if (abs(hi - lo) < 1e-9) {
    lo <- lo - 0.2
    hi <- hi + 0.2
  }
  seq(lo, hi, by = 0.2)
}

args <- commandArgs(trailingOnly = TRUE)
run_dir <- get_arg(args, "--run-dir", "")
sample_name <- get_arg(args, "--sample-name", "P20_LUAD")
score_csv <- get_arg(args, "--score-csv", "")
out_pdf <- get_arg(args, "--out-pdf", "")
out_png <- get_arg(args, "--out-png", "")
out_csv <- get_arg(args, "--out-csv", "")
root_dir <- get_arg(args, "--root-dir", getwd())

if (!nzchar(run_dir)) stop("--run-dir is required")
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
root_dir <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)

if (!nzchar(score_csv)) {
  score_csv <- file.path(root_dir, "results", "reproduce_paper", "tables", "spatial_scores", paste0(sample_name, "_scores.csv"))
}
if (!file.exists(score_csv)) stop("Missing score csv: ", score_csv)
score_csv <- normalizePath(score_csv, winslash = "/", mustWork = TRUE)

if (!nzchar(out_pdf)) {
  out_pdf <- file.path(run_dir, "figures", "Fig5B_spatial_highstem_activity_fromA.pdf")
}
if (!nzchar(out_png)) {
  out_png <- file.path(run_dir, "figures", "Fig5B_spatial_highstem_activity_fromA.png")
}
if (!nzchar(out_csv)) {
  out_csv <- file.path(run_dir, "tables", sprintf("spatial_highstem_activity_%s_fromA.csv", sample_name))
}

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

sample_dir <- find_sample_dir(root_dir, sample_name)
pos <- read_positions(sample_dir)
scores <- fread(score_csv)

if (all(c("Barcode", "Score", "X_pxl", "Y_pxl", "In_Tissue") %in% names(scores))) {
  scores2 <- data.table(
    barcode = as.character(scores$Barcode),
    Score = as.numeric(scores$Score),
    in_tissue = as.integer(scores$In_Tissue),
    imagecol = as.numeric(scores$X_pxl),
    imagerow = as.numeric(scores$Y_pxl)
  )
} else if (all(c("barcode", "Score", "in_tissue", "imagecol", "imagerow") %in% names(scores))) {
  scores2 <- data.table(
    barcode = as.character(scores$barcode),
    Score = as.numeric(scores$Score),
    in_tissue = as.integer(scores$in_tissue),
    imagecol = as.numeric(scores$imagecol),
    imagerow = as.numeric(scores$imagerow)
  )
} else {
  stop("Unsupported score csv columns in ", score_csv)
}

merged <- merge(
  scores2,
  pos[, .(barcode, in_tissue_pos = as.integer(in_tissue), row = as.integer(array_row), col = as.integer(array_col),
          pxl_row_in_fullres = as.numeric(pxl_row_in_fullres), pxl_col_in_fullres = as.numeric(pxl_col_in_fullres))],
  by = "barcode",
  all.x = TRUE
)

merged[, in_tissue := fifelse(!is.na(in_tissue), in_tissue, in_tissue_pos)]
merged[is.na(imagerow), imagerow := pxl_row_in_fullres]
merged[is.na(imagecol), imagecol := pxl_col_in_fullres]
merged <- merged[is.finite(Score) & is.finite(imagecol) & is.finite(imagerow)]
fwrite(
  merged[, .(barcode, Score, in_tissue, row, col, imagerow, imagecol)],
  out_csv
)

img_path <- file.path(sample_dir, "spatial", "tissue_lowres_image.png")
scale_path <- file.path(sample_dir, "spatial", "scalefactors_json.json")
if (!file.exists(img_path)) stop("Missing tissue_lowres_image.png in ", sample_dir)
if (!file.exists(scale_path)) stop("Missing scalefactors_json.json in ", sample_dir)

img <- readPNG(img_path)
img_faded <- make_rgba(img, alpha = 0.30)
sf <- fromJSON(scale_path)
h <- dim(img)[1]
w <- dim(img)[2]
tissue_sf <- as.numeric(sf$tissue_lowres_scalef)
spot_radius <- as.numeric(sf$spot_diameter_fullres) * tissue_sf / 2 * 0.92

plot_df <- copy(merged)
plot_df[, x := imagecol * tissue_sf]
plot_df[, y := h - imagerow * tissue_sf]
plot_df <- plot_df[in_tissue == 1]
hex_df <- hex_vertices_dt(plot_df, spot_radius)

score_breaks <- pretty_breaks_02(plot_df$Score)
score_limits <- range(score_breaks)
score_colors <- c("#2C7FB8", "#ABDDA4", "#FEE08B", "#F46D43", "#D73027")

img_grob <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
img_faded_grob <- rasterGrob(img_faded, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

p_left <- ggplot() +
  annotation_custom(img_grob, xmin = 0, xmax = w, ymin = 0, ymax = h) +
  coord_fixed(xlim = c(0, w), ylim = c(0, h), expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 4, 0, 0))

p_right <- ggplot() +
  annotation_custom(img_faded_grob, xmin = 0, xmax = w, ymin = 0, ymax = h) +
  geom_polygon(
    data = hex_df,
    aes(x = x, y = y, group = barcode, fill = Score),
    color = NA,
    alpha = 0.82,
    linewidth = 0
  ) +
  coord_fixed(xlim = c(0, w), ylim = c(0, h), expand = FALSE, clip = "off") +
  scale_fill_gradientn(
    colors = score_colors,
    limits = score_limits,
    breaks = score_breaks,
    oob = squish,
    name = "HighStem activity",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(22, "mm"),
      barheight = unit(6, "pt")
    )
  ) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8, face = "plain"),
    legend.text = element_text(size = 8),
    plot.margin = margin(0, 0, 0, 4)
  )

panel <- plot_grid(p_left, p_right, ncol = 2, rel_widths = c(1, 1))

pdf_w <- 90 / 25.4
pdf_h <- 45 / 25.4
ggsave(out_pdf, panel, width = pdf_w, height = pdf_h, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)

cat("OUT_PDF=", normalizePath(out_pdf, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("OUT_CSV=", normalizePath(out_csv, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("SAMPLE_DIR=", sample_dir, "\n", sep = "")
