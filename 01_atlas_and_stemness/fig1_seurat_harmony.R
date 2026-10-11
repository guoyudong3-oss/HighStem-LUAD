#!/usr/bin/env Rscript
scripts_needed <- c("Seurat", "ggplot2", "cowplot", "harmony", "patchwork", "KernSmooth", "sp", "scales")
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com", dependencies = TRUE)
  }
}
for (pkg in scripts_needed) install_if_missing(pkg)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(harmony)
  library(patchwork)
  nebulosa_available <- requireNamespace("Nebulosa", quietly = TRUE)
})

parse_args <- function(args) {
  run_dir <- NULL
  seurat_path <- NULL
  use_raw <- FALSE
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--run-dir") {
      run_dir <- args[i + 1]
      i <- i + 2
      next
    }
    if (args[i] == "--seurat") {
      seurat_path <- args[i + 1]
      i <- i + 2
      next
    }
    if (args[i] == "--raw") {
      use_raw <- TRUE
      i <- i + 1
      next
    }
    stop("Unknown argument: ", args[i])
  }
  if (is.null(run_dir)) {
    stop("Usage: Rscript fig1_seurat_harmony.R --run-dir <path> [--seurat <rds>] [--raw]")
  }
  list(
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    seurat_path = if (!is.null(seurat_path)) normalizePath(seurat_path, winslash = "/", mustWork = FALSE) else NULL,
    use_raw = use_raw
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
run_dir <- args$run_dir
set.seed(20260125)
fig_dir <- file.path(run_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

tab <- NULL
counts <- NULL

if (!is.null(args$seurat_path) && file.exists(args$seurat_path)) {
  cat("   Loading Seurat object: ", args$seurat_path, "\n")
  tab <- readRDS(args$seurat_path)
} else if (args$use_raw) {
  cat("   Loading raw UMI matrix from data/raw/GSE131907 ...\n")
  data_dir <- "data/raw/GSE131907"
  umi_matrix_file <- file.path(data_dir, "GSE131907_Lung_Cancer_raw_UMI_matrix.txt")
  annotation_file <- file.path(data_dir, "GSE131907_Lung_Cancer_cell_annotation.txt")
  if (!file.exists(umi_matrix_file) || !file.exists(annotation_file)) {
    stop("Missing raw GSE131907 files under data/raw/GSE131907/")
  }
  suppressPackageStartupMessages({
    library(data.table)
  })
  umi_data <- fread(umi_matrix_file, data.table = FALSE)
  gene_names <- umi_data[, 1]
  umi_data <- umi_data[, -1, drop = FALSE]
  rownames(umi_data) <- gene_names
  counts <- as(as.matrix(umi_data), "dgCMatrix")
  rm(umi_data); gc()
  meta_df <- read.table(annotation_file, header = TRUE, sep = "\t", row.names = 1)
  tab <- CreateSeuratObject(counts = counts, meta.data = meta_df, project = "PlasticStem_Lung")
} else {
  counts_path <- file.path("results/reproduce_paper/r_data", "full_matrix.mtx")
  features_path <- file.path("results/reproduce_paper/r_data", "full_features.tsv")
  barcodes_path <- file.path("results/reproduce_paper/r_data", "full_barcodes.tsv")
  meta_path <- file.path("results/reproduce_paper/r_data", "full_metadata.csv")

  missing_files <- c(counts_path, features_path, barcodes_path, meta_path)[
    !file.exists(c(counts_path, features_path, barcodes_path, meta_path))
  ]
  if (length(missing_files) > 0) {
    stop("Missing full dataset files: ", paste(missing_files, collapse = ", "))
  }

  read_mtx_auto <- function(path) {
    header <- readLines(path, n = 1)
    if (!grepl("^%%MatrixMarket", header)) {
      stop("Not a MatrixMarket file: ", path)
    }
    if (grepl(" array ", header)) {
      con <- file(path, open = "rt")
      on.exit(close(con))
      readLines(con, n = 1)
      line <- readLines(con, n = 1)
      while (startsWith(line, "%")) {
        line <- readLines(con, n = 1)
      }
      dims <- scan(text = line, what = integer(), quiet = TRUE)
      if (length(dims) < 2) stop("Invalid array MatrixMarket dims in ", path)
      n_row <- dims[1]
      n_col <- dims[2]
      vals <- scan(con, what = double(), n = n_row * n_col, quiet = TRUE)
      mat <- matrix(vals, nrow = n_row, ncol = n_col)
      return(mat)
    }
    Matrix::readMM(path)
  }

  cat("   Loading full counts matrix (all QC-passed cells)...\n")
  counts <- read_mtx_auto(counts_path)
  features <- read.table(features_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  gene_names <- make.unique(features[, 2])
  rownames(counts) <- gene_names
  barcodes <- readLines(barcodes_path)
  colnames(counts) <- barcodes

  cat("   Loading metadata...\n")
  meta_df <- read.csv(meta_path, stringsAsFactors = FALSE, row.names = 1)

  tab <- CreateSeuratObject(counts = counts, meta.data = meta_df, project = "PlasticStem_Lung")
}

marker_genes <- c(
  "NKG7", "KLRD1", "GNLY", "CD8A", "CD3D", "CD3E", "CCL5",
  "JCHAIN", "IGKC", "IGHG1",
  "EPCAM", "CDH1", "KRT7", "PGC", "PGA3", "MUC5AC",
  "PDGFRA", "LUM", "DCN",
  "CD79A", "CD79B", "MS4A1",
  "VWF", "PECAM1",
  "C1QA", "MARCO", "C1QB", "CD68", "LYZ", "FCN1", "S100A9", "S100A8",
  "RGS5", "PDGFRB", "ACTA2",
  "MKI67", "TOP2A", "CDK1",
  "CPA3", "GATA2", "KIT",
  "CD1C", "FCER1A", "CLEC9A",
  "CDH19", "PLP1", "SOX10"
)
rna_features <- rownames(tab[["RNA"]])
marker_genes <- marker_genes[marker_genes %in% rna_features]

# Nature (NPG) inspired palette
npg_base <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
  "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"
)

needs_processing <- args$use_raw
if (!needs_processing) {
  if (!("umap" %in% names(tab@reductions))) needs_processing <- TRUE
  if (!("harmony" %in% names(tab@reductions))) needs_processing <- TRUE
  if (!("seurat_clusters" %in% colnames(tab@meta.data))) needs_processing <- TRUE
}

if (needs_processing) {
  cat("Normalizing + PCA + Harmony...\n")
  min_val <- if (!is.null(counts)) suppressWarnings(min(counts)) else NA_real_
  if (is.finite(min_val) && min_val < 0) {
    warning("Detected negative values in full_matrix.mtx; skipping NormalizeData to avoid log1p NaNs.")
    counts_layer <- GetAssayData(tab, layer = "counts")
    tab <- SetAssayData(tab, layer = "data", new.data = counts_layer)
  } else {
    tab <- NormalizeData(tab)
  }
  tab <- tryCatch(
    FindVariableFeatures(tab, selection.method = "vst", nfeatures = 3000),
    error = function(e) {
      warning("FindVariableFeatures failed: ", conditionMessage(e))
      tab
    }
  )
  if (length(VariableFeatures(tab)) == 0) {
    warning("No variable features detected; using all features for scaling/PCA.")
    VariableFeatures(tab) <- rownames(tab)
  }
  tab <- ScaleData(tab, features = VariableFeatures(tab), verbose = FALSE)
  tab <- RunPCA(tab, features = VariableFeatures(tab), verbose = FALSE)
  tab <- RunHarmony(tab, group.by.vars = "Sample", assay.use = "RNA", plot_convergence = FALSE)
  tab <- RunUMAP(tab, reduction = "harmony", dims = 1:20)
  tab <- FindNeighbors(tab, reduction = "harmony", dims = 1:20)
  tab <- FindClusters(tab, resolution = 0.8)
} else {
  cat("Using existing Harmony/UMAP/clusters from Seurat object.\n")
}

cell_type_col <- "Cell_type.refined"
if (!cell_type_col %in% colnames(tab@meta.data)) {
  for (col in c("Cell_type", "Cell_subtype")) {
    if (col %in% colnames(tab@meta.data)) {
      cell_type_col <- col
      break
    }
  }
}

# Map subtypes to major cell classes to approximate the paper layout
subtypes <- if ("Cell_subtype" %in% colnames(tab@meta.data)) tab@meta.data$Cell_subtype else rep(NA, ncol(tab))
major <- rep(NA_character_, length(subtypes))
assign_major <- function(values, label) {
  major[subtypes %in% values] <<- label
}
assign_major(c("Naive CD4+ T", "CD4+ Th", "Treg", "CD8 low T", "Exhausted CD8+ T",
               "Naive CD8+ T", "Cytotoxic CD8+ T", "CD8+/CD4+ Mixed Th",
               "Exhausted Tfh"), "T")
assign_major(c("NK"), "NK")
assign_major(c("Follicular B cells", "MALT B cells", "GC B cells in the DZ",
               "GC B cells in the LZ", "GrB-secreting B cells"), "B")
assign_major(c("Plasma cells"), "Plasma")
assign_major(c("Monocytes"), "Mon")
assign_major(c("mo-Mac", "Alveolar Mac", "Pleural Mac", "Microglia/Mac"), "Mac")
assign_major(c("CD1c+ DCs", "CD163+CD14+ DCs", "Activated DCs", "CD141+ DCs",
               "CD207+CD1a+ LCs", "pDCs"), "MDC")
assign_major(c("MAST"), "Mast")
assign_major(c("Tip-like ECs", "Tumor ECs", "Stalk-like ECs", "Lymphatic ECs", "EPCs"), "Endo")
assign_major(c("Pericytes", "Smooth muscle cells"), "Mural")
assign_major(c("COL14A1+ matrix FBs", "COL13A1+ matrix FBs", "Myofibroblasts", "FB-like cells"), "Fib")
assign_major(c("AT2", "AT1", "Ciliated", "Club", "tS1", "tS2", "tS3",
               "Malignant cells", "Mesothelial cells"), "Epi")

fallback <- tab@meta.data[[cell_type_col]]
# GSE131907 has many NA values in Cell_type.refined / Cell_subtype; for major classes
# we fall back to the coarse Cell_type rather than treating them as "Other".
if ("Cell_type" %in% colnames(tab@meta.data)) {
  fallback[is.na(fallback)] <- tab@meta.data$Cell_type[is.na(fallback)]
}
fallback[is.na(fallback)] <- "Other"
major[is.na(major)] <- as.character(fallback[is.na(major)])
major[major == "T/NK cells"] <- "T"
major[major == "T lymphocytes"] <- "T"
major[major == "NK cells"] <- "NK"
major[major == "B lymphocytes"] <- "B"
major[major == "Fibroblasts"] <- "Fib"
major[major == "Endothelial cells"] <- "Endo"
major[major == "MAST cells"] <- "Mast"
major[major == "Epithelial cells"] <- "Epi"
major[major == "Myeloid cells"] <- "Mac"
major[major == "Oligodendrocytes"] <- "Glial"
major[major == "Undetermined"] <- "Other"
major <- factor(
  major,
  levels = c("NK", "T", "B", "Plasma", "MDC", "Mon", "Mac", "Mast", "Fib", "Mural", "Endo", "Epi", "Glial", "Other"),
  ordered = TRUE
)
major[is.na(major)] <- "Other"
tab@meta.data$celltype_major <- droplevels(major)

# Derive display labels (use dataset annotations; only refine NA/Other into Cyc/Glial)
celltype_display <- as.character(tab@meta.data$celltype_major)
celltype_display[celltype_display == "Other"] <- "Other"

cyc_genes <- c("MKI67", "TOP2A", "CDK1")
glial_genes <- c("PLP1", "SOX10", "CDH19")

# "Cyc" isn't provided by the dataset; approximate the paper by tagging cycling
# epithelial cells (so immune markers don't bleed into the Cyc column in DotPlot).
if (all(cyc_genes %in% rownames(tab))) {
  tab <- AddModuleScore(tab, features = list(cyc_genes), name = "CycScore")
  cyc_score <- tab@meta.data$CycScore1
  is_epi <- if ("Cell_type" %in% colnames(tab@meta.data)) {
    tab@meta.data$Cell_type == "Epithelial cells"
  } else {
    celltype_display == "Epi"
  }
  if (sum(is_epi, na.rm = TRUE) > 0) {
    thresh <- quantile(cyc_score[is_epi], 0.97, na.rm = TRUE)
    celltype_display[is_epi & cyc_score >= thresh] <- "Cyc"
  }
}

if (any(celltype_display == "Other") && all(glial_genes %in% rownames(tab))) {
  tab <- AddModuleScore(tab, features = list(glial_genes), name = "GlialScore")
  glial_score <- tab@meta.data$GlialScore1
  thresh <- quantile(glial_score[celltype_display == "Other"], 0.99, na.rm = TRUE)
  celltype_display[celltype_display == "Other" & glial_score >= thresh] <- "Glial"
}

tab@meta.data$celltype_display <- factor(
  celltype_display,
  levels = c("NK", "T", "B", "Plasma", "MDC", "Mon", "Mac", "Mast", "Fib", "Mural", "Endo", "Epi", "Cyc", "Glial", "Other"),
  ordered = TRUE
)

tab@meta.data$celltype_dot <- factor(
  celltype_display,
  levels = c("NK", "T", "Plasma", "Epi", "Fib", "B", "Endo", "Mac", "Mon", "Mural", "Cyc", "Mast", "MDC", "Glial"),
  ordered = TRUE
)

celltype_levels <- levels(tab@meta.data$celltype_display)
base_colors <- setNames(colorRampPalette(npg_base)(length(celltype_levels)), celltype_levels)
celltype_palette <- c(
  # Colors sampled from the original paper's Fig1C bottom strip (GSE131907 version).
  "NK" = "#4482B7",
  "T" = "#F17E24",
  "Plasma" = "#4AB035",
  "Epi" = "#DA2F2D",
  "Fib" = "#9367C6",
  "B" = "#8A554B",
  "Endo" = "#DE7ACE",
  "Mac" = "#B7C838",
  "Mon" = "#B6CC63",
  "Mural" = "#62D4D2",
  "Cyc" = "#B1CBE7",
  "Mast" = "#EEB778",
  "MDC" = "#A3F092",
  "Glial" = "#F39396",
  "Other" = "#CCCCCC"
)
celltype_colors <- base_colors
celltype_colors[names(celltype_palette)] <- celltype_palette[names(celltype_palette)]

cat("Saving panels...\n")
p_theme <- theme_minimal(base_size = 12) + theme(
  panel.grid = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9)
)

cluster_levels <- if ("seurat_clusters" %in% colnames(tab@meta.data)) {
  sort(unique(as.character(tab@meta.data$seurat_clusters)))
} else {
  character(0)
}
cluster_colors <- if (length(cluster_levels) > 0) {
  setNames(colorRampPalette(npg_base)(length(cluster_levels)), cluster_levels)
} else {
  NULL
}

pA <- DimPlot(
  tab,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  cols = cluster_colors
) +
  labs(title = "A. UMAP clusters") + p_theme

umap_coords <- Embeddings(tab, "umap")
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Fig1B (pragmatic fallback): Seurat/ggplot figure that conveys the same meaning
# (cell-type annotation + Normal/Tumor distribution) without reproducing the exact
# circular ring from the paper. Make it readable (bigger points, labels).
normal_origin <- c("nLung", "nLN")
tab@meta.data$origin_group <- ifelse(tab@meta.data$Sample_Origin %in% normal_origin, "Normal", "Tumor")
tab@meta.data$origin_group <- factor(tab@meta.data$origin_group, levels = c("Normal", "Tumor"), ordered = TRUE)

# With raster=TRUE, Seurat uses scattermore where pt.size is in raster pixels.
# If raster.dpi increases without scaling pt.size, points look too thin/faded.
b_raster_px <- 1000
b_pt_size <- 0.70 * (b_raster_px / 600)

pB_umap <- DimPlot(
  tab,
  reduction = "umap",
  group.by = "celltype_display",
  cols = celltype_colors,
  # Keep points fully opaque so colors match Fig1C palette.
  pt.size = b_pt_size,
  alpha = 1,
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  shuffle = TRUE,
  raster = TRUE,
  raster.dpi = c(b_raster_px, b_raster_px)
) +
  NoAxes() +
  labs(title = "B. Cell-type annotation") +
  guides(color = guide_legend(override.aes = list(size = 5.8, alpha = 1), ncol = 1)) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 10.2),
    legend.key.height = grid::unit(6.4, "mm"),
    legend.key.width = grid::unit(6.4, "mm"),
    legend.spacing.y = grid::unit(2.4, "mm"),
    legend.background = element_rect(fill = alpha("white", 0.92), color = NA),
    legend.position = c(0.95, 0.965),
    legend.justification = c("right", "top"),
    plot.margin = margin(5.5, 6, 5.5, 5.5)
  )

pB <- pB_umap

dot_cells <- !is.na(tab@meta.data$celltype_dot)
tab_dot <- tab[, dot_cells]
DefaultAssay(tab_dot) <- "RNA"
if (!"data" %in% SeuratObject::Layers(tab_dot[["RNA"]])) {
  tab_dot <- NormalizeData(tab_dot, verbose = FALSE)
}
Idents(tab_dot) <- tab_dot@meta.data$celltype_dot
dot_data <- DotPlot(tab_dot, features = marker_genes, group.by = "celltype_dot", assay = "RNA")$data

celltype_order <- levels(tab@meta.data$celltype_dot)
dot_data$id <- factor(dot_data$id, levels = celltype_order)
dot_data$features.plot <- factor(dot_data$features.plot, levels = rev(marker_genes))
dot_data$x <- as.numeric(dot_data$id)
dot_data$y <- as.numeric(dot_data$features.plot)
# Color: global scaling (not per-gene), with mild quantile stretching so the
# blue-white-red separation is closer to the paper.
q_from <- as.numeric(stats::quantile(dot_data$avg.exp.scaled, probs = c(0.02, 0.98), na.rm = TRUE))
dot_data$avg.exp.scaled.plot <- scales::rescale(dot_data$avg.exp.scaled, to = c(0, 1), from = q_from)
dot_data$avg.exp.scaled.plot <- pmin(pmax(dot_data$avg.exp.scaled.plot, 0), 1)
dot_data$avg.exp.scaled.plot <- dot_data$avg.exp.scaled.plot^0.7
dot_data$pct.exp <- pmin(pmax(dot_data$pct.exp, 0), 75)
dot_data$pct.exp_plot <- dot_data$pct.exp

dot_strip <- data.frame(
  x = seq_along(celltype_order),
  celltype = factor(celltype_order, levels = celltype_order)
)

gg_dot_main <- ggplot(dot_data, aes(x = x, y = y)) +
  geom_point(aes(size = pct.exp_plot, color = avg.exp.scaled.plot), na.rm = TRUE, stroke = 0) +
  scale_x_continuous(
    breaks = seq_along(celltype_order),
    labels = celltype_order,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq_along(marker_genes),
    labels = rev(marker_genes),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_color_gradientn(
    # Deeper maroon high end to match the paper's "red/purple" saturation.
    colors = c("#2166AC", "#F7F7F7", "#67001F"),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    name = "Relative\nExpression"
  ) +
  scale_size(
    # Larger max size so high-pct dots match the paper visually.
    range = c(0.10, 7.40),
    limits = c(0, 75),
    breaks = c(25, 50, 75),
    name = "Percentage\nExpressed"
  ) +
  coord_cartesian(
    # Add slight padding so large dots at the edges aren't clipped by the frame.
    xlim = c(0.3, length(celltype_order) + 0.7),
    ylim = c(0.3, length(marker_genes) + 0.7)
  ) +
  theme_bw(base_size = 11) +
  labs(title = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    # Paper-style bottom strip uses its own "lollipop" stems; avoid double-thick stems from axis ticks.
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length = grid::unit(0, "cm"),
    axis.text.y = element_text(face = "italic", size = 8.8),
    plot.title = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    # Give y-labels more room in combined figure to avoid left-edge clipping.
    plot.margin = margin(5.5, 14, 0, 34)
  ) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2, override.aes = list(color = "black"))
  )

make_dotplot_with_bottom_strip <- function(p, celltype_order, celltype_colors) {
  # Build a paper-style bottom strip using grid units (mm) so spacing is stable
  # and doesn't depend on ggplot coordinate scaling.
  stopifnot(all(celltype_order %in% names(celltype_colors)))
  n <- length(celltype_order)

  # Keep the strip compact but readable.
  strip_height <- grid::unit(22, "mm")
  # Match the paper: very short "sticks", circles nearly touching the labels.
  stem_len <- grid::unit(0.8, "mm")
  # Slightly smaller dots for a tighter, less crowded strip.
  dot_r <- grid::unit(2.6, "mm")
  label_gap <- grid::unit(0.1, "mm")

  xs <- (seq_len(n) - 0.5) / n

  y_top <- grid::unit(1, "npc")
  y_dot <- y_top - stem_len - dot_r
  y_label <- y_dot - dot_r - label_gap

  stems <- grid::segmentsGrob(
    x0 = grid::unit(xs, "npc"),
    x1 = grid::unit(xs, "npc"),
    y0 = y_top,
    y1 = y_dot + dot_r,
    gp = grid::gpar(col = "black", lwd = 0.7, lineend = "butt")
  )

  dots <- grid::circleGrob(
    x = grid::unit(xs, "npc"),
    y = y_dot,
    r = dot_r,
    gp = grid::gpar(fill = unname(celltype_colors[celltype_order]), col = NA)
  )

  labels <- grid::textGrob(
    label = celltype_order,
    x = grid::unit(xs, "npc"),
    y = y_label,
    rot = 90,
    # Anchor label top near the circles so the gap is minimal and consistent.
    just = c("right", "top"),
    gp = grid::gpar(col = "black", fontsize = 9.5, fontface = "bold")
  )

  strip <- grid::grobTree(stems, dots, labels)

  gt <- ggplot2::ggplotGrob(p)
  panel_idx <- which(gt$layout$name == "panel")
  panel_b <- gt$layout$b[panel_idx]
  panel_l <- gt$layout$l[panel_idx]
  panel_r <- gt$layout$r[panel_idx]

  gt <- gtable::gtable_add_rows(gt, heights = strip_height, pos = panel_b)
  gt <- gtable::gtable_add_grob(
    gt,
    strip,
    t = panel_b + 1,
    l = panel_l,
    r = panel_r,
    name = "bottom_strip",
    # Labels sit outside the panel frame; allow drawing beyond the gtable cell.
    clip = "off"
  )
  gt
}

gt_dot <- make_dotplot_with_bottom_strip(gg_dot_main, celltype_order, celltype_colors)

outline_idx <- chull(umap_coords[, 1], umap_coords[, 2])
outline_df <- data.frame(
  UMAP_1 = umap_coords[outline_idx, 1],
  UMAP_2 = umap_coords[outline_idx, 2]
)

density_genes <- c("NKG7", "JCHAIN", "EPCAM", "PDGFRA",
                   "CD79A", "VWF", "C1QA", "LYZ",
                   "RGS5", "MKI67", "CPA3", "CDH19")
density_genes <- density_genes[density_genes %in% rna_features]

make_density <- function(gene) {
  expr <- FetchData(tab, vars = gene, assay = "RNA", layer = "data")[, 1]
  df <- data.frame(umap_coords, expr = expr)
  df <- df[is.finite(df$expr) & df$expr > 0, , drop = FALSE]
  if (nebulosa_available) {
    p <- Nebulosa::plot_density(tab, genes = gene, reduction = "umap", scale = "log", pointsize = 0)
    p <- p + labs(title = gene) + theme_void()
  } else {
    if (nrow(df) < 50) {
      return(ggplot() + theme_void() + labs(title = gene))
    }
    if (nrow(df) > 60000) {
      df <- df[sample.int(nrow(df), 60000), , drop = FALSE]
    }
    rng_x <- range(umap_coords[, 1], na.rm = TRUE)
    rng_y <- range(umap_coords[, 2], na.rm = TRUE)
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
    # Mask outside the global UMAP outline to match the paper's black background.
    pip <- sp::point.in.polygon(grid$UMAP_1, grid$UMAP_2, outline_df$UMAP_1, outline_df$UMAP_2)
    grid$density[pip == 0] <- NA_real_
    # Rescale for a paper-like colorbar (0.00..0.20) regardless of UMAP coordinate range.
    max_d <- max(grid$density, na.rm = TRUE)
    if (is.finite(max_d) && max_d > 0) {
      grid$density <- grid$density / max_d * 0.20
    }
    brks <- pretty(c(0, max(grid$density, na.rm = TRUE)), n = 5)
    brks <- brks[brks >= 0]
    p <- ggplot(grid, aes(UMAP_1, UMAP_2, fill = density)) +
      geom_raster(interpolate = TRUE, na.rm = TRUE) +
      scale_fill_gradientn(
        colors = c("#0b0b0b", "#2d115f", "#5a1e78", "#93347a", "#d84d6a", "#f9915e", "#fcd67a"),
        na.value = "black",
        breaks = brks,
        labels = scales::number_format(accuracy = 0.01),
        name = "Density"
      ) +
      coord_equal(xlim = rng_x, ylim = rng_y) +
      theme_void() +
      labs(title = gene)
  }
  p <- p +
    geom_path(
      data = outline_df,
      aes(UMAP_1, UMAP_2),
      inherit.aes = FALSE,
      color = "white",
      linetype = "dashed",
      linewidth = 0.4
    ) +
    theme(
      plot.title = element_text(color = "white", size = 10, hjust = 0.5),
      legend.position = "right",
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      legend.key = element_rect(fill = "black", color = NA),
      legend.text = element_text(color = "white", size = 7),
      legend.title = element_text(color = "white", size = 7)
    )
  return(p)
}

density_plots <- lapply(density_genes, make_density)
pD <- wrap_plots(density_plots, ncol = 3)

pdf(file.path(fig_dir, "Fig1A_umap_clusters.pdf"), width = 6, height = 6)
print(pA)
dev.off()

pdf(file.path(fig_dir, "Fig1B_celltype_umap.pdf"), width = 6.2, height = 6.2)
print(pB)
dev.off()

pdf(file.path(fig_dir, "Fig1C_marker_dotplot.pdf"), width = 6.0, height = 10)
grid::grid.newpage()
grid::grid.draw(gt_dot)
dev.off()

pdf(file.path(fig_dir, "Fig1D_marker_density.pdf"), width = 12, height = 9)
print(pD)
dev.off()

# Publication-style combined Fig1 (A/B/C/D) with unified typography and compact layout.
pA_pub <- pA +
  labs(title = NULL) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(4, 4, 4, 4)
  )

pB_pub <- pB +
  labs(title = NULL) +
  theme(plot.margin = margin(4, 4, 4, 4))

gg_dot_pub <- gg_dot_main +
  labs(title = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8.0),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
gt_dot_pub <- make_dotplot_with_bottom_strip(gg_dot_pub, celltype_order, celltype_colors)
pC_pub <- patchwork::wrap_elements(full = gt_dot_pub)

pD_pub <- pD +
  plot_annotation(title = NULL) &
  theme(plot.margin = margin(0, 4, 4, 4))
# Add a thin white strip above D so the panel tag is always readable.
pD_pub <- (patchwork::plot_spacer() / pD_pub) + plot_layout(heights = c(0.055, 0.945))

fig1_core <- (
  # Keep top row close to the reference balance.
  pA_pub + pB_pub + plot_layout(widths = c(0.46, 0.54))
) / (
  # Constrained from reference + readability: C should not be over-squeezed.
  pC_pub + pD_pub + plot_layout(widths = c(0.34, 0.66))
) +
  # Reference-like row balance (bottom row dominant).
  plot_layout(heights = c(0.41, 0.59))

# Draw panel tags at figure level to avoid nested patchwork auto-tag artifacts.
fig1_pub <- cowplot::ggdraw(fig1_core) +
  cowplot::draw_plot_label(
    label = c("A", "B", "C", "D"),
    # Place tags at the start of each panel block according to measured splits.
    x = c(0.01, 0.465, 0.01, 0.345),
    y = c(0.995, 0.995, 0.605, 0.605),
    hjust = 0,
    vjust = 1,
    size = 18,
    fontface = "bold"
  )

ggsave(
  filename = file.path(fig_dir, "Fig1_combined_publication.pdf"),
  plot = fig1_pub,
  # Match reference figure aspect ratio (about 2000:2357).
  width = 10.2,
  height = 12.0,
  units = "in",
  limitsize = FALSE
)
ggsave(
  filename = file.path(fig_dir, "Fig1_combined_publication.png"),
  plot = fig1_pub,
  width = 10.2,
  height = 12.0,
  units = "in",
  dpi = 450,
  limitsize = FALSE
)

cat("Fig1 panels saved under ", fig_dir, "\n")
