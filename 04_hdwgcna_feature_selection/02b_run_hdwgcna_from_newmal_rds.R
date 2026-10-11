
library(Seurat)
library(hdWGCNA)
library(WGCNA)
library(Matrix)
library(data.table)
library(dplyr)

disableWGCNAThreads()

parse_args <- function(args) {
  run_dir <- NA_character_
  input_rds <- "results/reproduce_paper/runs/paper_main_epiref_full_noplot_20260209_114337/rds/seurat_malignant_with_stemness_runfinal_ref_mean_3sd.rds"
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--run-dir") {
      if (i == length(args)) {
        stop("Missing value for --run-dir")
      }
      run_dir <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--run-dir=")) {
      run_dir <- sub("^--run-dir=", "", a)
      i <- i + 1
      next
    }
    if (a == "--input-rds") {
      if (i == length(args)) {
        stop("Missing value for --input-rds")
      }
      input_rds <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--input-rds=")) {
      input_rds <- sub("^--input-rds=", "", a)
      i <- i + 1
      next
    }
    i <- i + 1
  }
  if (is.na(run_dir) || nchar(run_dir) == 0) {
    stop("Usage: Rscript 02b_run_hdwgcna_from_newmal_rds.R --run-dir <path> [--input-rds <seurat_rds>]")
  }
  list(
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    input_rds = normalizePath(input_rds, winslash = "/", mustWork = FALSE)
  )
}

cli <- parse_args(commandArgs(trailingOnly = TRUE))
run_dir <- cli$run_dir
input_rds <- cli$input_rds
message(paste("Using run dir:", run_dir))
message(paste("Using input rds:", input_rds))

dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(run_dir)) {
  stop(paste0("Failed to create run dir: ", run_dir))
}
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

tables_dir <- file.path(run_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

rds_dir <- file.path(run_dir, "rds")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

tom_dir <- file.path(run_dir, "tom")
dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

figures_dir <- file.path(run_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

out_hubs <- file.path(tables_dir, "wgcna_hub_genes.csv")
out_module_cor <- file.path(tables_dir, "wgcna_module_correlation.csv")
out_modules <- file.path(tables_dir, "wgcna_modules.csv")
out_power_table <- file.path(tables_dir, "wgcna_power_table.csv")
out_soft_power <- file.path(tables_dir, "wgcna_selected_soft_power.txt")
out_mes_cells <- file.path(tables_dir, "wgcna_module_eigengenes_cells.csv")
out_hubs_detail <- file.path(tables_dir, "wgcna_hub_genes_detail.csv")
out_hub_genes_all <- file.path(tables_dir, "wgcna_hub_genes_all.csv")

out_metacell_rds <- file.path(rds_dir, "HighStem_metacell_obj.rds")
out_metacell_meta <- file.path(tables_dir, "HighStem_metacell_metadata.csv")
out_cell_to_metacell <- file.path(tables_dir, "HighStem_cell_to_metacell.csv")

out_fig4a <- file.path(figures_dir, "Fig4A_soft_threshold.pdf")
out_fig4b <- file.path(figures_dir, "Fig4B_dendrogram.pdf")
out_fig4c <- file.path(figures_dir, "Fig4C_top_genes_by_module.pdf")
out_fig4d <- file.path(figures_dir, "Fig4D_module_umap.pdf")
out_fig4e <- file.path(figures_dir, "Fig4E_module_correlation.pdf")
out_fig4f <- file.path(figures_dir, "Fig4F_module_enrichment_dotplot.pdf")

if (!file.exists(input_rds)) {
  stop(paste0("Input rds not found: ", input_rds))
}

oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(tom_dir)

message("Loading malignant Seurat object from input rds ...")
seurat_obj <- readRDS(input_rds)
if (!inherits(seurat_obj, "Seurat")) {
  stop("input_rds is not a Seurat object: ", input_rds)
}
message(paste("Input object:", ncol(seurat_obj), "cells,", nrow(seurat_obj), "genes"))
if (!"cytotrace_score" %in% colnames(seurat_obj@meta.data)) {
  stop("Missing required metadata column cytotrace_score in input rds.")
}

group_col_candidates <- c("Sample", "sample_id")
group_col <- group_col_candidates[group_col_candidates %in% colnames(seurat_obj@meta.data)][1]
if (is.na(group_col) || nchar(group_col) == 0) {
  stop("Missing required sample grouping column in input rds. Tried: Sample, sample_id")
}

k_metacell <- 50
min_cells <- 51
target_metacells <- 1000
sample_counts <- table(seurat_obj@meta.data[[group_col]])
keep_samples <- names(sample_counts[sample_counts >= min_cells])
if (length(keep_samples) == 0) {
  stop(paste0("No groups in ", group_col, " have >= ", min_cells, " cells for metacell construction"))
}
if (length(keep_samples) < length(sample_counts)) {
  dropped <- setdiff(names(sample_counts), keep_samples)
  message(paste("Dropping small groups (<", min_cells, "cells):", paste(dropped, collapse = ", ")))
}
keep_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[group_col]] %in% keep_samples]
seurat_obj <- subset(seurat_obj, cells = keep_cells)
Idents(seurat_obj) <- group_col
message(paste("Processing malignant cells after group filter:", ncol(seurat_obj)))

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

message("Setting up hdWGCNA...")
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "HighStem"
)

message(paste0("Constructing Metacells by ", group_col, " (k=", k_metacell, ", min_cells=", min_cells, ", target_metacells=", target_metacells, ")..."))
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c(group_col),
  k = k_metacell,
  min_cells = min_cells,
  max_shared = 15,
  ident.group = group_col,
  target_metacells = target_metacells
)

seurat_obj <- NormalizeMetacells(seurat_obj)

metacell_obj <- GetMetacellObject(seurat_obj)
message(paste("DEBUG: Metacell Object Dimensions:", paste(dim(metacell_obj), collapse = " x ")))
if (ncol(metacell_obj) > 5000) {
  message("WARNING: Metacell count seems high. Check k value or grouping.")
}

if (!"cytotrace_score" %in% colnames(seurat_obj@meta.data)) {
  stop("Missing required column in seurat_obj@meta.data: cytotrace_score")
}

metacell_cols <- grep("metacell", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
message(paste("Metacell mapping columns in seurat_obj@meta.data:", paste(metacell_cols, collapse = ", ")))

metacell_col <- NULL
if ("metacell" %in% metacell_cols) {
  metacell_col <- "metacell"
} else if ("metacell_id" %in% metacell_cols) {
  metacell_col <- "metacell_id"
} else if (length(metacell_cols) > 0) {
  metacell_col <- metacell_cols[1]
}

if (is.null(metacell_col)) {
  stop("Could not find a metacell membership column (matching /metacell/i) in seurat_obj@meta.data")
}

cell_md <- seurat_obj@meta.data
cell_md$.__cell_id <- rownames(cell_md)
cell_md$.__metacell_id <- as.character(cell_md[[metacell_col]])

if (any(is.na(cell_md$.__metacell_id)) || any(cell_md$.__metacell_id == "")) {
  stop(paste("Metacell membership column has NA/empty values:", metacell_col))
}

write.csv(
  cell_md[, c(".__cell_id", ".__metacell_id", "cytotrace_score"), drop = FALSE],
  out_cell_to_metacell,
  row.names = FALSE
)
message(paste("Saved cell-to-metacell mapping to", out_cell_to_metacell))

metacell_cytotrace <- cell_md %>%
  group_by(.__metacell_id) %>%
  summarise(
    cytotrace_score = mean(cytotrace_score, na.rm = TRUE),
    stemness_score = mean(cytotrace_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

metacell_md <- metacell_obj@meta.data
metacell_md$.__metacell_id <- rownames(metacell_md)
metacell_md <- metacell_md %>% left_join(metacell_cytotrace, by = ".__metacell_id")

rownames(metacell_md) <- metacell_md$.__metacell_id
metacell_md$.__metacell_id <- NULL
metacell_obj@meta.data <- metacell_md

saveRDS(metacell_obj, out_metacell_rds)
write.csv(metacell_obj@meta.data, out_metacell_meta, row.names = TRUE)
message(paste("Saved metacell object to", out_metacell_rds))
message(paste("Saved metacell metadata to", out_metacell_meta))

if (!"cytotrace_score" %in% colnames(metacell_obj@meta.data)) {
  stop("Metacell metadata missing cytotrace_score after aggregation")
}

all_groups <- as.character(unique(Idents(metacell_obj)))
message(paste("Groups found in metacells:", paste(all_groups, collapse = ", ")))

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = all_groups,
  use_metacells = TRUE
)

message("Selecting Soft Power...")
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed"
)

power_table <- as.data.frame(GetPowerTable(seurat_obj))
write.csv(power_table, out_power_table, row.names = FALSE)
message(paste("Saved power table to", out_power_table))

valid_powers <- power_table[power_table$SFT.R.sq > 0.8, ]

if (nrow(valid_powers) > 0) {
  soft_power <- valid_powers$Power[1]
} else {
  message("Warning: No power reached 0.8 R2, using fallback 9")
  soft_power <- 9
}
message(paste("Selected Soft Power:", soft_power))
writeLines(as.character(soft_power), out_soft_power)
message(paste("Saved selected soft power to", out_soft_power))

message("Constructing Network (this takes time)...")
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = soft_power,
  setDatExpr = FALSE,
  tom_name = "HighStem_TOM",
  overwrite_tom = TRUE
)

seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

modules <- GetModules(seurat_obj)
write.csv(modules, out_modules, row.names = FALSE)
message(paste("Saved module assignments to", out_modules))

message("Correlating Modules with CytoTRACE (metacell-level MM/GS logic)...")
MEs <- GetMEs(seurat_obj, harmonized = TRUE)
write.csv(as.data.frame(MEs), out_mes_cells, row.names = TRUE)
message(paste("Saved module eigengenes (cell-level) to", out_mes_cells))

metacell_obj <- GetMetacellObject(seurat_obj)
if (!"cytotrace_score" %in% colnames(metacell_obj@meta.data)) {
  message("cytotrace_score missing in metacell metadata after network; rebuilding from cell mapping.")
  metacell_cytotrace_fix <- cell_md %>%
    group_by(.__metacell_id) %>%
    summarise(cytotrace_score = mean(cytotrace_score, na.rm = TRUE), .groups = "drop")

  metacell_md_fix <- metacell_obj@meta.data
  metacell_md_fix$.__metacell_id <- rownames(metacell_md_fix)
  metacell_md_fix <- metacell_md_fix %>% left_join(metacell_cytotrace_fix, by = ".__metacell_id")
  rownames(metacell_md_fix) <- metacell_md_fix$.__metacell_id
  metacell_md_fix$.__metacell_id <- NULL
  metacell_obj@meta.data <- metacell_md_fix
  seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj)
}

if (any(is.na(metacell_obj@meta.data$cytotrace_score)) && "cells_merged" %in% colnames(metacell_obj@meta.data)) {
  message("Found NA cytotrace_score in metacells; filling from cells_merged mapping.")
  sc_scores <- seurat_obj@meta.data$cytotrace_score
  names(sc_scores) <- rownames(seurat_obj@meta.data)
  cell_lists <- strsplit(as.character(metacell_obj@meta.data$cells_merged), ",", fixed = TRUE)
  filled_scores <- vapply(
    cell_lists,
    function(cell_ids) {
      vals <- sc_scores[cell_ids]
      if (all(is.na(vals))) {
        return(NA_real_)
      }
      mean(vals, na.rm = TRUE)
    },
    numeric(1)
  )
  na_before <- sum(is.na(metacell_obj@meta.data$cytotrace_score))
  metacell_obj@meta.data$cytotrace_score <- ifelse(
    is.na(metacell_obj@meta.data$cytotrace_score),
    filled_scores,
    metacell_obj@meta.data$cytotrace_score
  )
  na_after <- sum(is.na(metacell_obj@meta.data$cytotrace_score))
  message(paste("Metacell cytotrace_score NA:", na_before, "->", na_after))
  seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj)
}

trait_data <- metacell_obj@meta.data[["cytotrace_score"]]
if (is.null(trait_data)) {
  stop("cytotrace_score not found in metacell metadata.")
}
trait_vec <- as.numeric(trait_data)
names(trait_vec) <- rownames(metacell_obj@meta.data)

datExpr <- GetDatExpr(seurat_obj)
datExpr <- as.matrix(datExpr)
trait_vec <- trait_vec[rownames(datExpr)]
if (all(is.na(trait_vec))) {
  stop("Trait vector is all NA after aligning to datExpr rows (metacells).")
}
if (any(is.na(trait_vec))) {
  keep_mc <- !is.na(trait_vec)
  message(paste("Dropping", sum(!keep_mc), "metacells with NA cytotrace_score before correlation."))
  trait_vec <- trait_vec[keep_mc]
  datExpr <- datExpr[keep_mc, , drop = FALSE]
}

modules_nongrey <- modules %>% filter(module != "grey")
colors <- as.character(modules_nongrey$module[match(colnames(datExpr), modules_nongrey$gene_name)])
keep <- !is.na(colors)
datExpr_use <- datExpr[, keep, drop = FALSE]
colors_use <- colors[keep]
if (ncol(datExpr_use) == 0) {
  stop("No usable genes found in datExpr after module color mapping.")
}

MEs_mc <- WGCNA::moduleEigengenes(expr = datExpr_use, colors = colors_use)$eigengenes
MEs_mc <- as.data.frame(MEs_mc)
colnames(MEs_mc) <- gsub("^ME", "", colnames(MEs_mc))

cor_res <- cor(MEs_mc, trait_vec, use = "p")
p_res <- corPvalueStudent(cor_res, nrow(MEs_mc))

cor_df <- data.frame(
  module = rownames(cor_res),
  correlation = as.numeric(cor_res[, 1]),
  pval = as.numeric(p_res[, 1])
)
cor_df <- cor_df %>% arrange(desc(correlation))
print(head(cor_df))

top_module <- cor_df$module[1]
message(paste("Top Stemness Module:", top_module))

pos_modules <- cor_df %>% filter(correlation > 0.3 & pval < 0.05) %>% pull(module)
if (length(pos_modules) == 0) {
  pos_modules <- top_module
}

message(paste("Extracting genes from:", paste(pos_modules, collapse = ", ")))

mm_mat <- cor(datExpr_use, as.matrix(MEs_mc), use = "p")
mm_mat <- as.matrix(mm_mat)
rownames(mm_mat) <- colnames(datExpr_use)
mm_idx <- match(colors_use, colnames(mm_mat))
mm_vec <- rep(NA_real_, length(mm_idx))
valid_mm <- !is.na(mm_idx)
if (any(valid_mm)) {
  mm_vec[valid_mm] <- mm_mat[cbind(which(valid_mm), mm_idx[valid_mm])]
}

gs_vec <- cor(datExpr_use, trait_vec, use = "p")
gs_vec <- as.numeric(gs_vec)

hub_genes_selected <- data.frame(
  gene_name = colnames(datExpr_use),
  module = colors_use,
  MM = as.numeric(mm_vec),
  GS = gs_vec
) %>%
  filter(module %in% pos_modules) %>%
  filter(MM > 0.8 & GS > 0.2) %>%
  distinct()

hub_genes <- GetHubGenes(seurat_obj, n_hubs = 5000)
write.csv(hub_genes_selected, out_hubs_detail, row.names = FALSE)
message(paste("Saved selected hub genes detail to", out_hubs_detail))

selected_hubs <- hub_genes_selected %>% pull(gene_name)

write.csv(hub_genes, out_hub_genes_all, row.names = FALSE)
message(paste("Saved all hub genes to", out_hub_genes_all))

message(paste("Total Hub Genes selected:", length(selected_hubs)))

write.csv(data.frame(gene = selected_hubs), out_hubs, row.names = FALSE)
message(paste("Saved hub genes to", out_hubs))

write.csv(cor_df, out_module_cor, row.names = FALSE)
message(paste("Saved module correlations to", out_module_cor))

message("hdWGCNA analysis complete.")

message("Generating hdWGCNA figures...")

top_modules <- cor_df %>% arrange(desc(correlation)) %>% slice_head(n = 3) %>% pull(module)
if (length(top_modules) == 0) {
  top_modules <- pos_modules
}

preferred_modules <- c("brown", "green", "yellow")
modules_for_plot <- preferred_modules[preferred_modules %in% modules$module]
if (length(modules_for_plot) == 0) {
  modules_for_plot <- top_modules
}

res <- try({
  pdf(out_fig4a, width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(
    power_table$Power,
    power_table$SFT.R.sq,
    type = "b",
    xlab = "Soft Power",
    ylab = "Scale-free Topology Fit (R^2)",
    main = "Soft-thresholding"
  )
  abline(h = 0.8, col = "red", lty = 2)
  plot(
    power_table$Power,
    power_table$mean.k.,
    type = "b",
    xlab = "Soft Power",
    ylab = "Mean Connectivity",
    main = "Mean Connectivity"
  )
  dev.off()
}, silent = TRUE)
if (inherits(res, "try-error")) {
  message("Fig4A plot failed: ", res)
  if (dev.cur() != 1) dev.off()
}

res <- try({
  pdf(out_fig4b, width = 10, height = 6)
  PlotDendrogram(seurat_obj, main = "HighStem Dendrogram")
  dev.off()
}, silent = TRUE)
if (inherits(res, "try-error")) {
  message("Fig4B plot failed: ", res)
  if (dev.cur() != 1) dev.off()
}

res <- try({
  mm_col <- NULL
  if ("MM" %in% colnames(hub_genes)) {
    mm_col <- "MM"
  } else if ("kME" %in% colnames(hub_genes)) {
    mm_col <- "kME"
  }
  pdf(out_fig4c, width = 10, height = 7)
  for (m in top_modules) {
    df <- hub_genes[hub_genes$module == m, , drop = FALSE]
    if (nrow(df) == 0) next
    if (!is.null(mm_col)) {
      df <- df[order(df[[mm_col]], decreasing = TRUE), , drop = FALSE]
    }
    df <- head(df, 10)
    vals <- if (!is.null(mm_col)) df[[mm_col]] else seq_len(nrow(df))
    par(mar = c(8, 4, 3, 1))
    barplot(
      vals,
      names.arg = df$gene_name,
      las = 2,
      main = paste("Top genes in", m, "module"),
      ylab = ifelse(is.null(mm_col), "Rank", mm_col)
    )
  }
  dev.off()
}, silent = TRUE)
if (inherits(res, "try-error")) {
  message("Fig4C plot failed: ", res)
  if (dev.cur() != 1) dev.off()
}

res <- try({
  me_cols <- c()
  for (m in top_modules) {
    if (paste0("ME", m) %in% colnames(MEs)) {
      me_cols <- c(me_cols, paste0("ME", m))
    } else if (m %in% colnames(MEs)) {
      me_cols <- c(me_cols, m)
    }
  }
  if (length(me_cols) > 0) {
    common_cells <- intersect(rownames(seurat_obj@meta.data), rownames(MEs))
    seurat_obj@meta.data[common_cells, me_cols] <- MEs[common_cells, me_cols]
    pdf(out_fig4d, width = 4 * length(me_cols), height = 4)
    print(FeaturePlot(seurat_obj, features = me_cols, reduction = "umap", ncol = length(me_cols)))
    dev.off()
  } else {
    message("Fig4D plot skipped: no matching module eigengene columns.")
  }
}, silent = TRUE)
if (inherits(res, "try-error")) {
  message("Fig4D plot failed: ", res)
  if (dev.cur() != 1) dev.off()
}

res <- try({
  me_mat <- cor(MEs, use = "p")
  text_mat <- signif(me_mat, 2)
  pdf(out_fig4e, width = 8, height = 8)
  labeledHeatmap(
    Matrix = me_mat,
    xLabels = colnames(me_mat),
    yLabels = colnames(me_mat),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = text_mat,
    setStdMargins = FALSE,
    cex.text = 0.6
  )
  dev.off()
}, silent = TRUE)
if (inherits(res, "try-error")) {
  message("Fig4E plot failed: ", res)
  if (dev.cur() != 1) dev.off()
}

res <- try({
  if ("cytotrace_score" %in% colnames(seurat_obj@meta.data)) {
    ct <- seurat_obj@meta.data$cytotrace_score
    q1 <- as.numeric(quantile(ct, 0.25, na.rm = TRUE))
    q3 <- as.numeric(quantile(ct, 0.75, na.rm = TRUE))
    stem_group <- ifelse(ct >= q3, "HighStem", ifelse(ct <= q1, "LowStem", "DTStem"))
    seurat_obj$stem_group_quartile <- stem_group

    module_genes_list <- lapply(modules_for_plot, function(m) {
      modules$gene_name[modules$module == m]
    })
    valid_idx <- lengths(module_genes_list) > 0
    module_genes_list <- module_genes_list[valid_idx]
    module_names <- modules_for_plot[valid_idx]
    if (length(module_genes_list) > 0) {
      seurat_obj <- AddModuleScore(seurat_obj, features = module_genes_list, name = "ModuleScore_")
      score_cols <- paste0("ModuleScore_", seq_along(module_genes_list))
      for (i in seq_along(score_cols)) {
        new_name <- paste0("ModuleScore_", module_names[i])
        seurat_obj@meta.data[[new_name]] <- seurat_obj@meta.data[[score_cols[i]]]
      }

      pdf(out_fig4f, width = 8, height = 5)
      print(
        DotPlot(
          seurat_obj,
          features = paste0("ModuleScore_", module_names),
          group.by = "stem_group_quartile"
        ) + RotatedAxis()
      )
      dev.off()
    } else {
      message("Fig4F plot skipped: no module genes available for plotting.")
    }
  } else {
    message("Fig4F plot skipped: cytotrace_score not found.")
  }
}, silent = TRUE)
if (inherits(res, "try-error")) {
  message("Fig4F plot failed: ", res)
  if (dev.cur() != 1) dev.off()
}

message("Figure generation complete.")
