#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(KernSmooth)
  library(sp)
  library(scales)
  library(grid)
})

parse_args <- function(args) {
  opt <- list(
    seurat = "results/reproduce_paper/r_data/seurat_obj.rds",
    out_dir = "results/reproduce_paper/runs/significant_noclobber_full/figures/fig1_manual_relayout_20260218",
    prefix = "Fig1D_marker_density_rebuilt",
    width_in = 7.0866,
    height_in = 7.9,
    dpi = 450
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key == "--seurat") {
      opt$seurat <- args[i + 1]; i <- i + 2; next
    }
    if (key == "--out-dir") {
      opt$out_dir <- args[i + 1]; i <- i + 2; next
    }
    if (key == "--prefix") {
      opt$prefix <- args[i + 1]; i <- i + 2; next
    }
    if (key == "--width-in") {
      opt$width_in <- as.numeric(args[i + 1]); i <- i + 2; next
    }
    if (key == "--height-in") {
      opt$height_in <- as.numeric(args[i + 1]); i <- i + 2; next
    }
    if (key == "--dpi") {
      opt$dpi <- as.integer(args[i + 1]); i <- i + 2; next
    }
    stop("Unknown arg: ", key)
  }
  opt
}

ensure_umap <- function(tab) {
  need <- !("umap" %in% names(tab@reductions))
  if (!need) return(tab)

  message("UMAP missing; recomputing Harmony + UMAP for panel D ...")
  if (!"data" %in% SeuratObject::Layers(tab[["RNA"]])) {
    tab <- NormalizeData(tab, verbose = FALSE)
  }
  tab <- FindVariableFeatures(tab, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  tab <- ScaleData(tab, features = VariableFeatures(tab), verbose = FALSE)
  tab <- RunPCA(tab, features = VariableFeatures(tab), verbose = FALSE)
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("harmony package required to recompute panel D UMAP.")
  }
  group_var <- if ("Sample" %in% colnames(tab@meta.data)) "Sample" else colnames(tab@meta.data)[1]
  tab <- harmony::RunHarmony(tab, group.by.vars = group_var, assay.use = "RNA", plot_convergence = FALSE)
  tab <- RunUMAP(tab, reduction = "harmony", dims = 1:20, verbose = FALSE)
  tab
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)

tab <- readRDS(opts$seurat)
tab <- ensure_umap(tab)

umap_coords <- Embeddings(tab, "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
rng_x <- range(umap_coords[, 1], na.rm = TRUE)
rng_y <- range(umap_coords[, 2], na.rm = TRUE)

outline_idx <- chull(umap_coords[, 1], umap_coords[, 2])
outline_df <- data.frame(
  UMAP_1 = umap_coords[outline_idx, 1],
  UMAP_2 = umap_coords[outline_idx, 2]
)

genes <- c("NKG7", "JCHAIN", "EPCAM",
           "PDGFRA", "CD79A", "VWF",
           "C1QA", "LYZ", "RGS5",
           "MKI67", "CPA3", "CDH19")
genes <- genes[genes %in% rownames(tab)]

make_density <- function(gene) {
  expr <- FetchData(tab, vars = gene, assay = "RNA", layer = "data")[, 1]
  df <- data.frame(UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2], expr = expr)
  df <- df[is.finite(df$expr) & df$expr > 0, , drop = FALSE]

  if (nrow(df) < 50) {
    p <- ggplot() + theme_void()
    p <- p + theme(
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
    return(p)
  }
  if (nrow(df) > 60000) {
    set.seed(20260218)
    df <- df[sample.int(nrow(df), 60000), , drop = FALSE]
  }

  bw_x <- diff(rng_x) / 60
  bw_y <- diff(rng_y) / 60
  dens <- KernSmooth::bkde2D(
    as.matrix(df[, c("UMAP_1", "UMAP_2")]),
    bandwidth = c(bw_x, bw_y),
    gridsize = c(220, 220),
    range.x = list(rng_x, rng_y)
  )
  grid <- expand.grid(UMAP_1 = dens$x1, UMAP_2 = dens$x2)
  grid$density <- as.vector(dens$fhat)
  pip <- sp::point.in.polygon(grid$UMAP_1, grid$UMAP_2, outline_df$UMAP_1, outline_df$UMAP_2)
  grid$density[pip == 0] <- NA_real_

  max_d <- max(grid$density, na.rm = TRUE)
  if (!is.finite(max_d) || max_d <= 0) max_d <- 1e-6
  brks <- pretty(c(0, max_d), n = 3)
  brks <- brks[brks >= 0]
  lbl_acc <- 0.01  # copy: every colour bar uses 2 decimals

  p <- ggplot(grid, aes(UMAP_1, UMAP_2, fill = density)) +
    geom_raster(interpolate = TRUE, na.rm = TRUE) +
    scale_fill_gradientn(
      colors = c("#0b0b0b", "#2d115f", "#5a1e78", "#93347a", "#d84d6a", "#f9915e", "#fcd67a"),
      na.value = "black",
      limits = c(0, max_d),
      breaks = brks,
      labels = scales::number_format(accuracy = lbl_acc),
      name = "Density"
    ) +
    coord_equal(xlim = rng_x, ylim = rng_y, expand = FALSE) +
    theme_void() +
    geom_path(
      data = outline_df,
      aes(UMAP_1, UMAP_2),
      inherit.aes = FALSE,
      color = "white",
      linetype = "dashed",
      linewidth = 0.45
    ) +
    labs(title = gene) +   # copy: gene name above the black panel (not on top of the dashed outline)
    theme(
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA), legend.key.width = unit(2.5, "mm"), legend.key.height = unit(5.5, "mm"),
      legend.title = element_text(color = "black", size = 8),
      legend.text = element_text(color = "black", size = 8),
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5, color = "black", family = "Arial"),
      plot.margin = margin(0, 8, 6, 6)
    )
  p
}

plist <- lapply(genes, make_density)
pD <- wrap_plots(plist, ncol = 3, byrow = TRUE) +
  plot_layout(guides = "keep") +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  )

out_pdf <- file.path(opts$out_dir, paste0(opts$prefix, ".pdf"))
out_png <- file.path(opts$out_dir, paste0(opts$prefix, ".png"))

ggsave(out_pdf, pD, width = opts$width_in, height = opts$height_in, units = "in", dpi = opts$dpi, bg = "white", limitsize = FALSE, device = cairo_pdf)
ggsave(out_png, pD, width = opts$width_in, height = opts$height_in, units = "in", dpi = opts$dpi, bg = "white", limitsize = FALSE)

writeLines(c(out_pdf, out_png))
