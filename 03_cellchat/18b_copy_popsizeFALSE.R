#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(KEGGREST)
})

parse_args <- function(args) {
  get_arg <- function(key, default = NULL) {
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
  has_flag <- function(key) {
    if (any(args == key)) return(TRUE)
    key_eq <- paste0("^", key, "=")
    kv <- grep(key_eq, args, value = TRUE)
    if (length(kv) == 0) return(FALSE)
    v <- tolower(sub(key_eq, "", kv[length(kv)]))
    v %in% c("1", "true", "yes", "y")
  }
  list(
    run_dir = get_arg("--run-dir", NULL),
    malignant_rds = get_arg(
      "--malignant-rds",
      "results/reproduce_paper/runs/paper_main_epiref_full_noplot_20260209_114337/rds/seurat_malignant_with_stemness_runfinal_ref_mean_3sd.rds"
    ),
    skip_hi = has_flag("--skip-hi")
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$run_dir) || !nzchar(args$run_dir)) {
  stop("Usage: Rscript 18b_fig3_cellchat_pathway_metabolism.R --run-dir <path> [--malignant-rds <path>]")
}

run_dir <- normalizePath(args$run_dir, winslash = "/", mustWork = FALSE)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

fig_dir <- file.path(run_dir, "figures")
tab_dir <- file.path(run_dir, "tables")
ckpt_dir <- file.path(run_dir, "checkpoints_fig3")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

input_dir <- file.path(run_dir, "cellchat_input")
expr_file <- file.path(input_dir, "cellchat_expr_gene_by_cell.csv.gz")
meta_file <- file.path(input_dir, "cellchat_metadata_sampled.csv")

if (!file.exists(expr_file)) stop("Missing CellChat expr input: ", expr_file)
if (!file.exists(meta_file)) stop("Missing CellChat meta input: ", meta_file)
if (!isTRUE(args$skip_hi) && !file.exists(args$malignant_rds)) {
  stop("Missing malignant rds: ", args$malignant_rds)
}

resume_mode <- tolower(Sys.getenv("FIG3_RESUME", unset = "1")) %in% c("1", "true", "yes")

step_file <- function(id, name) {
  file.path(ckpt_dir, sprintf("%02d_%s.rds", as.integer(id), name))
}

restore_or_run <- function(id, name, expr) {
  p <- step_file(id, name)
  if (resume_mode && file.exists(p)) {
    cat("↪ Resume:", p, "\n")
    return(readRDS(p))
  }
  cat("▶ Step", sprintf("%02d", as.integer(id)), name, "\n")
  obj <- force(expr)
  saveRDS(obj, p)
  cat("✔ Saved:", p, "\n")
  obj
}

save_plot_both <- function(plot_obj, stem, width = 7, height = 5, dpi = 320) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot_obj, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot_obj, width = width, height = height, dpi = dpi)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_data_layer <- function(seu, assay = "RNA") {
  d1 <- tryCatch(SeuratObject::GetAssayData(seu, layer = "data", assay = assay), error = function(e) NULL)
  if (!is.null(d1) && nrow(d1) > 0 && ncol(d1) > 0) return(d1)
  seu <- Seurat::NormalizeData(seu, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  d2 <- tryCatch(SeuratObject::GetAssayData(seu, layer = "data", assay = assay), error = function(e) NULL)
  if (!is.null(d2) && nrow(d2) > 0 && ncol(d2) > 0) return(d2)
  stop("Cannot access RNA data layer from Seurat object.")
}

cat("================================================================================\n")
cat("🧭 Fig3 method-flow analysis (CellChat + Pathway + Metabolism)\n")
cat("================================================================================\n")
cat("run_dir:", run_dir, "\n")
cat("resume_mode:", ifelse(resume_mode, "ON", "OFF"), "\n")
cat("skip_hi:", ifelse(isTRUE(args$skip_hi), "ON", "OFF"), "\n")

# 01 Load CellChat inputs
inp <- restore_or_run(1, "load_cellchat_inputs", {
  expr_df <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_file, row.names = 1, check.names = FALSE)
  common <- intersect(colnames(expr_df), rownames(meta))
  expr_df <- expr_df[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]
  meta$cellchat_label <- as.character(meta$cellchat_label)
  list(expr = as.matrix(expr_df), meta = meta)
})

# 02 Create CellChat object
cellchat <- restore_or_run(2, "cellchat_create", {
  cc <- createCellChat(object = inp$expr, meta = inp$meta, group.by = "cellchat_label")
  cc@DB <- CellChatDB.human
  cc
})

# 03 subsetData
cellchat <- restore_or_run(3, "cellchat_subset_db", {
  subsetData(cellchat)
})

# 04 identify OE
cellchat <- restore_or_run(4, "cellchat_identify_overexpressed", {
  cc <- identifyOverExpressedGenes(cellchat)
  cc <- identifyOverExpressedInteractions(cc)
  cc
})

# 05 communication probability
cellchat <- restore_or_run(5, "cellchat_compute_prob", {
  cc <- computeCommunProb(cellchat, population.size = FALSE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc
})

# 06 pathway aggregation
cellchat <- restore_or_run(6, "cellchat_compute_pathway", {
  cc <- computeCommunProbPathway(cellchat)
  cc <- aggregateNet(cc)
  cc
})

# 06b centrality (required for signaling-role heatmaps/scatter)
if (!isTRUE(args$skip_hi)) {
  cellchat <- restore_or_run(61, "cellchat_compute_centrality", {
    netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  })
} else {
  cat("ℹ skip_hi=ON: skipping centrality and signaling-role plots (Fig3B/C).\n")
}

# 07 export base tables
cellchat_tables <- restore_or_run(7, "cellchat_export_tables", {
  df_net <- subsetCommunication(cellchat)
  write.csv(df_net, file.path(tab_dir, "fig3_cellchat_interactions_all.csv"), row.names = FALSE)
  count_mat <- cellchat@net$count
  weight_mat <- cellchat@net$weight
  write.csv(as.data.frame(count_mat), file.path(tab_dir, "fig3_cellchat_count_matrix.csv"), row.names = TRUE)
  write.csv(as.data.frame(weight_mat), file.path(tab_dir, "fig3_cellchat_weight_matrix.csv"), row.names = TRUE)
  list(
    n_interactions = nrow(df_net),
    labels = levels(cellchat@idents)
  )
})

# Shared color map across Fig3A/B/C so node/bar/scatter colors are consistent.
get_fig3_group_colors <- function(groups) {
  base_colors <- c(
    "B lymphocytes" = "#E41A1C",
    "DTStem" = "#377EB8",
    "Endothelial cells" = "#4DAF4A",
    "Fibroblasts" = "#984EA3",
    "HighStem" = "#F18F01",
    "LowStem" = "#F781BF",
    "MAST cells" = "#B39BC8",
    "Myeloid cells" = "#A65628",
    "NK cells" = "#56B4E9",
    "NormalEpi" = "#1F2A7C",
    "T lymphocytes" = "#1B9E77",
    "T/NK cells" = "#9BCF6E",
    "Undetermined" = "#D9B300"
  )
  miss <- setdiff(groups, names(base_colors))
  if (length(miss) > 0) {
    extra <- grDevices::hcl.colors(length(miss), palette = "Dynamic")
    names(extra) <- miss
    base_colors <- c(base_colors, extra)
  }
  base_colors[groups]
}

group_levels <- levels(cellchat@idents)
group_colors <- get_fig3_group_colors(group_levels)

# 08 Fig3A: circles (count + weight)
restore_or_run(8, "fig3A_cellchat_circle", {
  group_size <- as.numeric(table(cellchat@idents))
  weight_label_mult <- suppressWarnings(as.numeric(Sys.getenv("FIG3A_WEIGHT_LABEL_MULT", unset = "1000")))
  if (!is.finite(weight_label_mult) || weight_label_mult <= 0) weight_label_mult <- 1000
  w_plot <- cellchat@net$weight * weight_label_mult
  w_title <- if (abs(weight_label_mult - 1) > 1e-12) {
    paste0("Interaction weights/strength (labels x", format(weight_label_mult, scientific = FALSE), ")")
  } else {
    "Interaction weights/strength"
  }
  pdf(file.path(fig_dir, "Fig3A_cellchat_circle.pdf"), width = 14, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(
    cellchat@net$count,
    color.use = group_colors,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.label.cex = 0.65,
    title.name = "Number of interactions"
  )
  netVisual_circle(
    w_plot,
    color.use = group_colors,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.label.cex = 0.65,
    title.name = w_title
  )
  dev.off()
  png(file.path(fig_dir, "Fig3A_cellchat_circle.png"), width = 4200, height = 1800, res = 300)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(
    cellchat@net$count,
    color.use = group_colors,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.label.cex = 0.65,
    title.name = "Number of interactions"
  )
  netVisual_circle(
    w_plot,
    color.use = group_colors,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = TRUE,
    edge.label.cex = 0.65,
    title.name = w_title
  )
  dev.off()
  TRUE
})

if (!isTRUE(args$skip_hi)) {
  # 09 Fig3B: signaling role heatmaps (outgoing/incoming)
  restore_or_run(9, "fig3B_signaling_role_heatmaps", {
    fig3b_font_size <- suppressWarnings(as.numeric(Sys.getenv("FIG3B_FONT_SIZE", unset = "3.8")))
    if (!is.finite(fig3b_font_size) || fig3b_font_size <= 0) fig3b_font_size <- 3.8
    fig3b_font_title <- suppressWarnings(as.numeric(Sys.getenv("FIG3B_TITLE_FONT_SIZE", unset = "10")))
    if (!is.finite(fig3b_font_title) || fig3b_font_title <= 0) fig3b_font_title <- 10
    fig3b_top_pathways <- suppressWarnings(as.integer(Sys.getenv("FIG3B_TOP_PATHWAYS", unset = "40")))
    if (!is.finite(fig3b_top_pathways) || fig3b_top_pathways <= 0) fig3b_top_pathways <- 40L
    fig3b_ht_width <- suppressWarnings(as.numeric(Sys.getenv("FIG3B_HEATMAP_WIDTH", unset = "12")))
    if (!is.finite(fig3b_ht_width) || fig3b_ht_width <= 0) fig3b_ht_width <- 12
    fig3b_ht_height <- suppressWarnings(as.numeric(Sys.getenv("FIG3B_HEATMAP_HEIGHT", unset = "12")))
    if (!is.finite(fig3b_ht_height) || fig3b_ht_height <= 0) fig3b_ht_height <- 12

    signaling_use <- NULL
    prob_arr <- cellchat@netP$prob
    if (!is.null(prob_arr) && length(dim(prob_arr)) == 3) {
      pathway_strength <- apply(prob_arr, 3, sum, na.rm = TRUE)
      pathway_strength <- sort(pathway_strength, decreasing = TRUE)
      if (length(pathway_strength) > fig3b_top_pathways) {
        signaling_use <- names(pathway_strength)[seq_len(fig3b_top_pathways)]
      } else {
        signaling_use <- names(pathway_strength)
      }
    }

    ht_out <- netAnalysis_signalingRole_heatmap(
      cellchat,
      signaling = signaling_use,
      pattern = "outgoing",
      color.use = group_colors,
      width = fig3b_ht_width,
      height = fig3b_ht_height,
      font.size = fig3b_font_size,
      font.size.title = fig3b_font_title
    )
    ht_in <- netAnalysis_signalingRole_heatmap(
      cellchat,
      signaling = signaling_use,
      pattern = "incoming",
      color.use = group_colors,
      width = fig3b_ht_width,
      height = fig3b_ht_height,
      font.size = fig3b_font_size,
      font.size.title = fig3b_font_title
    )
    pdf(file.path(fig_dir, "Fig3B_signaling_role_heatmaps.pdf"), width = 14, height = 6)
    draw(ht_out + ht_in, ht_gap = grid::unit(0.5, "cm"))
    dev.off()
    png(file.path(fig_dir, "Fig3B_signaling_role_heatmaps.png"), width = 4200, height = 1800, res = 300)
    draw(ht_out + ht_in, ht_gap = grid::unit(0.5, "cm"))
    dev.off()
    TRUE
  })

  # 10 Fig3C: outgoing vs incoming scatter
  restore_or_run(10, "fig3C_signaling_role_scatter", {
    wmat <- as.matrix(cellchat@net$weight)
    cmat <- as.matrix(cellchat@net$count)
    diag(wmat) <- 0
    diag(cmat) <- 0
    d <- data.frame(
      labels = rownames(wmat),
      x = rowSums(wmat),
      y = colSums(wmat),
      Count = rowSums(cmat) + colSums(cmat),
      stringsAsFactors = FALSE
    ) %>% filter(is.finite(x), is.finite(y), x > 0, y > 0)
    d$labels <- factor(d$labels, levels = group_levels)

    # Clip only extreme tails to avoid a single dominant point compressing the panel.
    qx <- as.numeric(quantile(d$x, probs = 0.95, na.rm = TRUE))
    qy <- as.numeric(quantile(d$y, probs = 0.95, na.rm = TRUE))
    d <- d %>%
      mutate(
        x_plot = pmin(x, qx),
        y_plot = pmin(y, qy)
      )

    p <- ggplot(d, aes(x = x_plot, y = y_plot, color = labels, size = Count, label = labels)) +
      geom_point(alpha = 0.9) +
      ggrepel::geom_text_repel(size = 3.3, max.overlaps = 100, show.legend = FALSE) +
      scale_color_manual(values = group_colors, drop = FALSE) +
      scale_size_continuous(range = c(2.5, 11)) +
      labs(
        x = "Outgoing interaction strength",
        y = "Incoming interaction strength"
      ) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 11)
      )
    save_plot_both(p, "Fig3C_outgoing_incoming_scatter", 7, 6)
    TRUE
  })
}

short_partner_label <- function(x) {
  map <- c(
    "B lymphocytes" = "B",
    "Endothelial cells" = "Endo",
    "Fibroblasts" = "Fib",
    "Myeloid cells" = "Myeloid",
    "MAST cells" = "Mast",
    "NK cells" = "NK",
    "NormalEpi" = "Normalepi",
    "T lymphocytes" = "T",
    "T/NK cells" = "T/NK",
    "Undetermined" = "Undetermined",
    "DTStem" = "DTStem",
    "LowStem" = "LowStem",
    "HighStem" = "HighStem"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

build_compare_df <- function(
  count_mat,
  high = "HighStem",
  ref = "LowStem",
  mode = c("outgoing", "incoming"),
  exclude_groups = c("HighStem", "LowStem", "DTStem")
) {
  mode <- match.arg(mode)
  if (!(high %in% rownames(count_mat) && ref %in% rownames(count_mat))) return(data.frame())

  partners <- intersect(rownames(count_mat), colnames(count_mat))
  partners <- setdiff(partners, exclude_groups)

  df <- lapply(partners, function(p) {
    if (mode == "outgoing") {
      high_count <- as.numeric(count_mat[high, p])   # partner provides receptor
      ref_count  <- as.numeric(count_mat[ref, p])
    } else {
      high_count <- as.numeric(count_mat[p, high])   # partner provides ligand
      ref_count  <- as.numeric(count_mat[p, ref])
    }
    data.frame(
      partner = p,
      high_count = high_count,
      ref_count = ref_count,
      log2_ratio = log2((high_count + 1) / (ref_count + 1)),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  df <- df %>% filter((high_count + ref_count) > 0)
  if (nrow(df) == 0) return(df)

  df$partner_show <- short_partner_label(df$partner)
  # Order to mimic paper style: larger HighStem interaction counts closer to the bottom.
  ord <- df %>% arrange(desc(high_count), desc(ref_count), partner_show)
  df$partner_show <- factor(df$partner_show, levels = rev(unique(ord$partner_show)))

  ratio_clip <- suppressWarnings(as.numeric(Sys.getenv("FIG3DE_RATIO_CLIP", unset = "0.6")))
  if (!is.finite(ratio_clip) || ratio_clip <= 0) ratio_clip <- 0.6
  df$log2_ratio_clip <- pmax(pmin(df$log2_ratio, ratio_clip), -ratio_clip)
  df
}

plot_compare_bar <- function(df, mode = c("outgoing", "incoming")) {
  mode <- match.arg(mode)
  if (nrow(df) == 0) return(NULL)

  # Visual style target: paper-like mirrored bars + top log2 axis.
  # Keep the style defaults but expand to cover observed ratios so red-dot positions remain faithful.
  ratio_lims_base <- if (mode == "outgoing") c(-0.35, 0.55) else c(-0.05, 0.45)
  ratio_show <- df$log2_ratio
  rmin <- suppressWarnings(min(ratio_show, na.rm = TRUE))
  rmax <- suppressWarnings(max(ratio_show, na.rm = TRUE))
  rspan <- rmax - rmin
  if (!is.finite(rspan) || rspan <= 0) rspan <- max(1, abs(rmax), abs(rmin), na.rm = TRUE)
  rpad <- max(0.06 * rspan, 0.05)
  ratio_lims <- c(
    min(ratio_lims_base[1], rmin - rpad, 0),
    max(ratio_lims_base[2], rmax + rpad, 0)
  )
  if (!all(is.finite(ratio_lims)) || ratio_lims[2] <= ratio_lims[1]) ratio_lims <- ratio_lims_base

  x_left <- -max(df$ref_count, na.rm = TRUE) * 1.12
  x_right <- max(df$high_count, na.rm = TRUE) * 1.10
  if (!is.finite(x_left) || x_left >= 0) x_left <- -1
  if (!is.finite(x_right) || x_right <= 0) x_right <- 1

  # Linear mapping: ratio_lims[1] -> x_left, 0 -> 0, ratio_lims[2] -> x_right (piecewise so both limits are exact).
  ratio_to_x <- function(r) {
    ifelse(r <= 0, r * x_left / ratio_lims[1], r * x_right / ratio_lims[2])
  }
  ratio_clip <- pmin(pmax(ratio_show, ratio_lims[1]), ratio_lims[2])
  df$dot_x <- ratio_to_x(ratio_clip)
  x_pad <- 0.02 * (x_right - x_left)
  df$dot_x <- pmax(df$dot_x, x_left + x_pad)
  df$dot_x <- pmin(df$dot_x, x_right - x_pad)
  txt_inset_max <- 0.018 * (x_right - x_left)
  df$ref_label_x <- -df$ref_count + pmin(pmax(df$ref_count * 0.16, 0.18), txt_inset_max)
  df$high_label_x <- df$high_count - pmin(pmax(df$high_count * 0.16, 0.18), txt_inset_max)
  # Dot color strictly follows direction in data: red if HighStem > LowStem, else black.
  df$dot_col <- ifelse(df$high_count > df$ref_count, "red3", "black")

  y_label <- if (mode == "outgoing") {
    "Cell types interacting with Malignant\n(providing receptor)"
  } else {
    "Cell types interacting with Malignant\n(providing ligand)"
  }
  panel_letter <- if (mode == "outgoing") "D" else "E"

  # Bottom x-axis labels keep the paper style: hide negative-side labels.
  bot_breaks <- pretty(c(x_left, x_right), n = 5)
  bot_labels <- ifelse(bot_breaks < 0,
                       "",
                       as.character(as.integer(round(bot_breaks))))

  # Top secondary axis: log2 ratio ticks mapped to x-axis native coordinates.
  top_breaks_ratio <- pretty(ratio_lims, n = 4)
  top_breaks_ratio <- sort(unique(c(top_breaks_ratio, 0)))
  top_breaks_ratio <- top_breaks_ratio[
    top_breaks_ratio >= ratio_lims[1] - 1e-12 & top_breaks_ratio <= ratio_lims[2] + 1e-12
  ]
  if (length(top_breaks_ratio) < 2) {
    top_breaks_ratio <- sort(unique(c(ratio_lims[1], 0, ratio_lims[2])))
  }
  top_breaks_x     <- ratio_to_x(top_breaks_ratio)
  top_labels_ratio <- sapply(top_breaks_ratio, function(v) {
    if (abs(v) < 1e-10) return("0")
    sub("\\.?0+$", "", formatC(v, format = "f", digits = 2))
  })

  top_axis_title <- "Log2(# of interaction pairs in HighStem/\n# of interaction pairs in LowStem)"

  p <- ggplot(df, aes(y = partner_show)) +
    geom_col(aes(x = -ref_count), fill = "#38807F", color = "#38807F", linewidth = 0.9, width = 0.7) +
    geom_col(aes(x = high_count), fill = "#AB2329", color = "#AB2329", linewidth = 0.9, width = 0.7) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.5) +
    geom_text(
      data = df %>% filter(ref_count > 0),
      aes(x = ref_label_x, label = as.integer(ref_count)),
      color = "white", size = 4.5, fontface = "bold"
    ) +
    geom_text(
      data = df %>% filter(high_count > 0),
      aes(x = high_label_x, label = as.integer(high_count)),
      color = "white", size = 4.8, fontface = "bold"
    ) +
    geom_point(aes(x = dot_x, color = dot_col), size = 2.9, alpha = 0.95, show.legend = FALSE) +
    scale_color_identity() +
    scale_x_continuous(
      name = "# of interaction pairs",
      limits = c(x_left, x_right),
      breaks = bot_breaks,
      labels = bot_labels,
      sec.axis = sec_axis(
        ~ .,
        name   = top_axis_title,
        breaks = top_breaks_x,
        labels = top_labels_ratio
      )
    ) +
    scale_y_discrete(position = "right") +
    labs(y = y_label, tag = panel_letter) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.9),
      axis.line = element_blank(),
      axis.title.x.bottom = element_text(size = 11, family = "sans"),
      axis.title.x.top    = element_text(size = 10, family = "sans"),
      axis.text.x.bottom  = element_text(size = 10, family = "sans"),
      axis.text.x.top     = element_text(size = 9,  family = "sans"),
      axis.text.y         = element_text(size = 13, family = "sans"),
      axis.title.y        = element_text(size = 13, family = "sans"),
      plot.margin         = margin(18, 14, 8, 20),
      plot.tag            = element_text(family = "serif", face = "bold", size = 30),
      plot.tag.position   = c(0.01, 0.98)
    )
  p
}

plot_lr_bubble_custom <- function(df, out_stem, top_lr = 70L, w = 10.5, h = 7) {
  if (nrow(df) == 0) return(FALSE)

  # Keep significant interactions first; fallback to all if too sparse.
  df_sig <- df %>% filter(!is.na(pval), pval < 0.05)
  if (nrow(df_sig) >= 20) {
    df <- df_sig
  }

  lr_rank <- df %>%
    group_by(interaction_name_2) %>%
    summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_prob))
  keep_lr <- head(lr_rank$interaction_name_2, n = min(top_lr, nrow(lr_rank)))
  df <- df %>% filter(interaction_name_2 %in% keep_lr)
  if (nrow(df) == 0) return(FALSE)

  top_points <- suppressWarnings(as.integer(Sys.getenv("FIG3FG_TOP_POINTS", unset = "180")))
  if (is.na(top_points) || top_points < 50) top_points <- 180L
  df <- df %>% arrange(desc(prob)) %>% head(top_points)

  df <- df %>%
    mutate(
      pair = paste0(source, "->", target),
      pval_safe = pmax(pval, 1e-300),
      pscore = -log10(pval_safe)
    )

  pcap <- suppressWarnings(as.numeric(quantile(df$pscore, probs = 0.95, na.rm = TRUE)))
  pfloor <- suppressWarnings(as.numeric(quantile(df$pscore, probs = 0.05, na.rm = TRUE)))
  if (!is.finite(pcap)) pcap <- max(df$pscore, na.rm = TRUE)
  if (!is.finite(pfloor)) pfloor <- min(df$pscore, na.rm = TRUE)
  if (pcap <= pfloor) pcap <- pfloor + 1
  df$pscore_plot <- (pmin(pmax(df$pscore, pfloor), pcap) - pfloor) / (pcap - pfloor + 1e-12)

  lprob <- log10(pmax(df$prob, 1e-300))
  lmin <- min(lprob, na.rm = TRUE)
  lmax <- max(lprob, na.rm = TRUE)
  if (is.finite(lmin) && is.finite(lmax) && lmax > lmin) {
    df$prob_plot <- (lprob - lmin) / (lmax - lmin)
  } else {
    df$prob_plot <- 0.5
  }

  pair_order <- df %>%
    group_by(pair) %>%
    summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_prob)) %>%
    pull(pair)
  lr_order <- df %>%
    group_by(interaction_name_2) %>%
    summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_prob)) %>%
    pull(interaction_name_2)

  df$pair <- factor(df$pair, levels = pair_order)
  df$interaction_name_2 <- factor(df$interaction_name_2, levels = rev(lr_order))

  p <- ggplot(df, aes(x = pair, y = interaction_name_2, color = prob_plot, size = pscore_plot)) +
    geom_point(alpha = 0.9) +
    scale_color_gradientn(colors = c("#2C7BB6", "#66C2A5", "#FEE08B", "#F46D43", "#D73027")) +
    scale_size_continuous(range = c(1.4, 4.6), name = "-log10(p)") +
    labs(x = NULL, y = NULL, color = "Commun. Prob.") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_line(size = 0.2, color = "#e6e6e6"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  save_plot_both(p, out_stem, w, h)
  TRUE
}

plot_lr_bubble_paperstyle <- function(df, out_stem, mode = c("outgoing", "incoming"), panel_letter = NULL, w = 9.0, h = 10.0) {
  mode <- match.arg(mode)
  if (nrow(df) == 0) return(FALSE)

  p_cut <- suppressWarnings(as.numeric(Sys.getenv("FIG3FG_P_CUT", unset = "0.01")))
  if (!is.finite(p_cut) || p_cut <= 0 || p_cut >= 1) p_cut <- 0.01

  df_sig <- df %>% filter(!is.na(pval), pval <= p_cut)
  sig_label <- paste0("p <", format(p_cut, trim = TRUE))
  if (nrow(df_sig) < 30) {
    df_sig <- df %>% filter(!is.na(pval), pval <= 0.05)
    sig_label <- "p <0.05"
  }
  if (nrow(df_sig) == 0) {
    df_sig <- df
    sig_label <- "all"
  }

  top_lr <- suppressWarnings(as.integer(Sys.getenv("FIG3FG_TOP_LR_PAPER", unset = "48")))
  top_pairs <- suppressWarnings(as.integer(Sys.getenv("FIG3FG_TOP_PAIRS_PAPER", unset = "48")))
  if (!is.finite(top_lr) || top_lr < 10) top_lr <- 48L
  if (!is.finite(top_pairs) || top_pairs < 10) top_pairs <- 48L

  df_sig <- df_sig %>% mutate(pair = paste0(source, " -> ", target))
  keep_lr <- df_sig %>%
    group_by(interaction_name_2) %>%
    summarise(score = max(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score)) %>%
    slice_head(n = top_lr) %>%
    pull(interaction_name_2)
  keep_pair <- df_sig %>%
    group_by(pair) %>%
    summarise(score = max(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score)) %>%
    slice_head(n = top_pairs) %>%
    pull(pair)

  d <- df_sig %>%
    filter(interaction_name_2 %in% keep_lr, pair %in% keep_pair)
  if (nrow(d) == 0) return(FALSE)

  mal <- c("DTStem", "HighStem", "LowStem")
  if (mode == "outgoing") {
    src_ord <- mal[mal %in% unique(d$source)]
    tgt_ord <- sort(unique(d$target))
    ord_df <- expand.grid(source = src_ord, target = tgt_ord, stringsAsFactors = FALSE)
  } else {
    src_ord <- sort(unique(d$source))
    tgt_ord <- mal[mal %in% unique(d$target)]
    ord_df <- expand.grid(source = src_ord, target = tgt_ord, stringsAsFactors = FALSE)
  }
  ord_df$pair <- paste0(ord_df$source, " -> ", ord_df$target)
  pair_levels <- ord_df$pair[ord_df$pair %in% unique(d$pair)]
  pair_levels <- unique(c(pair_levels, setdiff(keep_pair, pair_levels), setdiff(unique(d$pair), pair_levels)))
  d$pair <- factor(d$pair, levels = pair_levels)

  lr_levels <- d %>%
    group_by(interaction_name_2) %>%
    summarise(score = max(prob, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score)) %>%
    pull(interaction_name_2)
  d$interaction_name_2 <- factor(d$interaction_name_2, levels = rev(lr_levels))

  pminv <- min(d$prob, na.rm = TRUE)
  pmaxv <- max(d$prob, na.rm = TRUE)
  if (!is.finite(pminv) || !is.finite(pmaxv) || pmaxv <= pminv) {
    d$prob_plot <- 0.5
  } else {
    d$prob_plot <- (d$prob - pminv) / (pmaxv - pminv)
  }
  d$sig_class <- sig_label

  p <- ggplot(d, aes(x = pair, y = interaction_name_2, color = prob_plot, size = sig_class)) +
    geom_point(alpha = 0.95) +
    scale_color_gradientn(
      colors = c("#3B4CC0", "#2C7FB8", "#41B6C4", "#A1D99B", "#FEE08B", "#F46D43", "#D73027"),
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("min", "max"),
      name = "Commun. Prob."
    ) +
    scale_size_manual(values = setNames(2.2, sig_label), name = "p-value") +
    labs(x = NULL, y = NULL, tag = panel_letter %||% "") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_line(size = 0.22, color = "#e8e8e8"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.8),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8.5),
      axis.text.y = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.tag = element_text(family = "serif", face = "bold", size = 34),
      plot.tag.position = c(0.01, 0.99)
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2, override.aes = list(color = "black"))
    )

  save_plot_both(p, out_stem, width = w, height = h)
  TRUE
}

# 11 Fig3D: HighStem vs LowStem (partner providing receptor; malignant outgoing)
restore_or_run(11, "fig3D_high_vs_low", {
  df <- build_compare_df(cellchat@net$count, high = "HighStem", ref = "LowStem", mode = "outgoing")
  write.csv(df, file.path(tab_dir, "fig3D_high_vs_low_interaction_pairs.csv"), row.names = FALSE)
  p <- plot_compare_bar(df, mode = "outgoing")
  if (!is.null(p)) save_plot_both(p, "Fig3D_highstem_vs_lowstem_pairs", 9.2, 7.0)
  TRUE
})

# 12 Fig3E: HighStem vs LowStem (partner providing ligand; malignant incoming)
restore_or_run(12, "fig3E_high_vs_low_incoming", {
  df <- build_compare_df(cellchat@net$count, high = "HighStem", ref = "LowStem", mode = "incoming")
  write.csv(df, file.path(tab_dir, "fig3E_high_vs_low_incoming_interaction_pairs.csv"), row.names = FALSE)
  # keep legacy filename for backward compatibility
  write.csv(df, file.path(tab_dir, "fig3E_high_vs_dtstem_interaction_pairs.csv"), row.names = FALSE)
  p <- plot_compare_bar(df, mode = "incoming")
  if (!is.null(p)) {
    save_plot_both(p, "Fig3E_highstem_vs_lowstem_pairs", 9.2, 7.0)
    save_plot_both(p, "Fig3E_highstem_vs_dtstem_pairs", 9.2, 7.0)
  }
  TRUE
})

# 13 Fig3F: bubble (HighStem -> others)
restore_or_run(13, "fig3F_bubble_outgoing", {
  df_net <- subsetCommunication(cellchat)
  write.csv(df_net, file.path(tab_dir, "fig3F_source_table_all.csv"), row.names = FALSE)
  df <- df_net %>%
    filter(source == "HighStem", target != "HighStem")
  write.csv(df, file.path(tab_dir, "fig3F_source_table_highstem_outgoing.csv"), row.names = FALSE)
  mal <- c("DTStem", "HighStem", "LowStem")
  df_plot <- df_net %>%
    filter(source %in% mal, !(target %in% mal))
  write.csv(df_plot, file.path(tab_dir, "fig3F_source_table_malignant_outgoing.csv"), row.names = FALSE)
  plot_lr_bubble_paperstyle(df_plot, out_stem = "Fig3F_bubble_highstem_outgoing", mode = "outgoing", panel_letter = "F", w = 9.2, h = 10.2)
  TRUE
})

# 14 Fig3G: bubble (others -> HighStem)
restore_or_run(14, "fig3G_bubble_incoming", {
  df_net <- subsetCommunication(cellchat)
  df <- df_net %>%
    filter(target == "HighStem", source != "HighStem")
  write.csv(df, file.path(tab_dir, "fig3G_source_table_highstem_incoming.csv"), row.names = FALSE)
  mal <- c("DTStem", "HighStem", "LowStem")
  df_plot <- df_net %>%
    filter(!(source %in% mal), target %in% mal)
  write.csv(df_plot, file.path(tab_dir, "fig3G_source_table_malignant_incoming.csv"), row.names = FALSE)
  plot_lr_bubble_paperstyle(df_plot, out_stem = "Fig3G_bubble_highstem_incoming", mode = "incoming", panel_letter = "G", w = 9.2, h = 10.2)
  TRUE
})

if (!isTRUE(args$skip_hi)) {
  # 15 load malignant object for H/I
  mal_data <- restore_or_run(15, "load_malignant_for_HI", {
    seu <- readRDS(args$malignant_rds)
    if (!"stem_group" %in% colnames(seu@meta.data)) stop("malignant rds missing stem_group.")
    expr <- get_data_layer(seu, assay = "RNA")
    list(seu = seu, expr = expr)
  })

# Helpers for KEGG gene sets
kegg_symbol_cache <- new.env(parent = emptyenv())
kegg_name_cache <- new.env(parent = emptyenv())
local_kegg_map <- NULL
kegg_pathway_name_map <- NULL
if (requireNamespace("scPagwas", quietly = TRUE)) {
  ok <- try({
    data("Genes_by_pathway_kegg", package = "scPagwas", envir = environment())
    if (exists("Genes_by_pathway_kegg", inherits = FALSE)) {
      local_kegg_map <- lapply(Genes_by_pathway_kegg, function(x) unique(toupper(as.character(x))))
    }
  }, silent = TRUE)
  if (inherits(ok, "try-error")) local_kegg_map <- NULL
}

load_kegg_pathway_name_map <- function(cache_file = file.path(tab_dir, "kegg_pathway_name_map.tsv")) {
  parse_lines <- function(lines) {
    if (is.null(lines) || length(lines) == 0) return(NULL)
    sp <- strsplit(lines, "\t", fixed = TRUE)
    ok <- vapply(sp, function(x) length(x) >= 2, logical(1))
    sp <- sp[ok]
    if (length(sp) == 0) return(NULL)
    ids <- vapply(sp, function(x) x[[1]], character(1))
    nms <- vapply(sp, function(x) x[[2]], character(1))
    nms <- sub(" - Homo sapiens \\(human\\)$", "", nms)
    out <- nms
    names(out) <- ids
    out
  }

  if (file.exists(cache_file)) {
    lines <- tryCatch(readLines(cache_file, warn = FALSE), error = function(e) NULL)
    mp <- parse_lines(lines)
    if (!is.null(mp) && length(mp) > 0) return(mp)
  }

  lines <- tryCatch(readLines("https://rest.kegg.jp/list/pathway/hsa", warn = FALSE), error = function(e) NULL)
  mp <- parse_lines(lines)
  if (!is.null(mp) && length(mp) > 0) {
    try(writeLines(paste(names(mp), mp, sep = "\t"), con = cache_file), silent = TRUE)
    return(mp)
  }
  NULL
}
kegg_pathway_name_map <- load_kegg_pathway_name_map()

safe_kegg_get <- function(path_id, retries = 4L) {
  for (i in seq_len(retries)) {
    obj <- tryCatch(keggGet(path_id), error = function(e) e)
    if (!inherits(obj, "error") && !is.null(obj) && length(obj) > 0) {
      return(obj[[1]])
    }
    Sys.sleep(min(2 * i, 6))
  }
  NULL
}
get_kegg_symbols <- function(path_id) {
  if (exists(path_id, envir = kegg_symbol_cache, inherits = FALSE)) {
    return(get(path_id, envir = kegg_symbol_cache, inherits = FALSE))
  }
  symbols <- character(0)
  if (!is.null(local_kegg_map) && path_id %in% names(local_kegg_map)) {
    symbols <- local_kegg_map[[path_id]]
  } else {
    obj1 <- safe_kegg_get(path_id)
    if (!is.null(obj1) && !is.null(obj1$GENE)) {
      gvec <- obj1$GENE
      if (length(gvec) >= 2) {
        symbols <- unique(toupper(sub(";.*$", "", gvec[seq(2, length(gvec), by = 2)])))
      }
    } else {
      message("KEGG fetch failed for: ", path_id)
    }
  }
  assign(path_id, symbols, envir = kegg_symbol_cache)
  symbols
}

get_kegg_name <- function(path_id) {
  if (!is.null(kegg_pathway_name_map) && path_id %in% names(kegg_pathway_name_map)) {
    return(kegg_pathway_name_map[[path_id]])
  }
  if (exists(path_id, envir = kegg_name_cache, inherits = FALSE)) {
    return(get(path_id, envir = kegg_name_cache, inherits = FALSE))
  }
  nm <- path_id
  obj <- safe_kegg_get(path_id)
  if (!is.null(obj) && !is.null(obj$NAME) && length(obj$NAME) > 0) {
    nm <- as.character(obj$NAME[[1]])
    nm <- sub(" - Homo sapiens \\(human\\)$", "", nm)
    nm <- sub(" \\- reference pathway$", "", nm)
  }
  assign(path_id, nm, envir = kegg_name_cache)
  nm
}

score_pathway_by_cells <- function(expr_mat, pathway_symbols) {
  rns <- toupper(rownames(expr_mat))
  idx <- which(rns %in% pathway_symbols)
  if (length(idx) < 5) {
    return(NULL)
  }
  colMeans(expr_mat[idx, , drop = FALSE])
}

residualize_vector <- function(score_vec, cov_vec, min_n = 30L) {
  ok <- is.finite(score_vec) & is.finite(cov_vec)
  if (sum(ok) < min_n) return(score_vec)
  if (stats::sd(cov_vec[ok]) < 1e-8) return(score_vec)
  fit <- stats::lm(score_vec[ok] ~ cov_vec[ok])
  out <- score_vec
  out[ok] <- stats::residuals(fit)
  out
}

  # 16 Fig3H: pathway activity heatmap (oncogenic signaling)
  restore_or_run(16, "fig3H_pathway_activity", {
  expr <- mal_data$expr
  meta <- mal_data$seu@meta.data
  groups <- as.character(meta$stem_group)
  names(groups) <- rownames(meta)
  common_cells <- intersect(colnames(expr), names(groups))
  expr <- expr[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  groups <- groups[common_cells]

  use_ncount_resid <- tolower(Sys.getenv("FIG3H_RESIDUALIZE_NCOUNT", unset = "1")) %in% c("1", "true", "yes")
  use_androgen_fallback <- tolower(Sys.getenv("FIG3H_ANDROGEN_FALLBACK", unset = "1")) %in% c("1", "true", "yes")
  ncount_cov <- NULL
  if (use_ncount_resid && "nCount_RNA" %in% colnames(meta)) {
    ncount_cov <- log1p(as.numeric(meta$nCount_RNA))
  }
  androgen_fallback_syms <- c(
    "AR", "NCOA1", "NCOA2", "NCOA3", "NCOR1", "NCOR2", "FOXA1", "FKBP5",
    "TMPRSS2", "KLK2", "KLK3", "PPARGC1A", "PIK3CA", "PIK3CB", "AKT1",
    "MAPK1", "MAPK3", "EGFR", "SRC", "PTEN", "CTNNB1", "GSK3B", "FOXO1"
  )

  oncogenic_kegg <- c(
    "PI3K" = "hsa04151",
    "WNT" = "hsa04310",
    "TGFb" = "hsa04350",
    "JAK-STAT" = "hsa04630",
    "EGFR" = "hsa04012",
    "MAPK" = "hsa04010",
    "TNFa" = "hsa04668",
    "VEGF" = "hsa04370",
    "p53" = "hsa04115",
    "Androgen" = "hsa04980",
    "Estrogen" = "hsa04915",
    "Hypoxia" = "hsa04066",
    "Trail" = "hsa04210",
    "NFkB" = "hsa04064"
  )

  score_list <- list()
  gene_count <- integer(0)
  debug_rows <- list()
  for (nm in names(oncogenic_kegg)) {
    syms <- get_kegg_symbols(oncogenic_kegg[[nm]])
    used_fallback <- FALSE
    if (nm == "Androgen" && use_androgen_fallback) {
      overlap0 <- length(intersect(toupper(rownames(expr)), syms))
      if (length(syms) == 0 || overlap0 < 5) {
        syms <- unique(toupper(c(syms, androgen_fallback_syms)))
        used_fallback <- TRUE
      }
    }
    overlap <- length(intersect(toupper(rownames(expr)), syms))
    debug_rows[[length(debug_rows) + 1]] <- data.frame(
      pathway = nm,
      kegg_id = oncogenic_kegg[[nm]],
      n_symbols = length(syms),
      n_overlap = overlap,
      used_fallback = used_fallback,
      stringsAsFactors = FALSE
    )
    scr <- score_pathway_by_cells(expr, syms)
    if (!is.null(scr)) {
      if (!is.null(ncount_cov)) {
        scr <- residualize_vector(scr, ncount_cov, min_n = 30L)
      }
      score_list[[nm]] <- scr
      gene_count[[nm]] <- overlap
    }
  }
  debug_df <- bind_rows(debug_rows)
  write.csv(debug_df, file.path(tab_dir, "fig3H_kegg_overlap_debug.csv"), row.names = FALSE)
  if (length(score_list) == 0) {
    stop("Fig3H: no valid pathways scored. See fig3H_kegg_overlap_debug.csv")
  }
  score_mat <- do.call(cbind, score_list)
  score_df <- as.data.frame(score_mat)
  score_df$stem_group <- groups
  agg <- score_df %>%
    group_by(stem_group) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop")

  group_order <- c("DTStem", "LowStem", "HighStem")
  keep_groups <- group_order[group_order %in% agg$stem_group]
  agg <- agg %>%
    mutate(stem_group = factor(stem_group, levels = keep_groups)) %>%
    arrange(stem_group)
  agg_mat <- t(as.matrix(agg[, setdiff(colnames(agg), "stem_group"), drop = FALSE]))
  colnames(agg_mat) <- as.character(agg$stem_group)

  z_mat <- t(scale(t(agg_mat)))
  z_mat[is.na(z_mat)] <- 0

  z_clip_q <- suppressWarnings(as.numeric(Sys.getenv("FIG3H_Z_CLIP_Q", unset = "0.9")))
  if (!is.finite(z_clip_q) || z_clip_q <= 0.5 || z_clip_q >= 1) z_clip_q <- 0.9
  z_abs <- suppressWarnings(as.numeric(stats::quantile(abs(as.vector(z_mat)), probs = z_clip_q, na.rm = TRUE)))
  if (!is.finite(z_abs) || z_abs <= 0) z_abs <- 0.6
  z_abs <- max(0.35, z_abs)
  z_draw <- pmax(pmin(z_mat, z_abs), -z_abs)

  write.csv(agg_mat, file.path(tab_dir, "fig3H_pathway_activity_group_mean.csv"), row.names = TRUE)
  write.csv(z_mat, file.path(tab_dir, "fig3H_pathway_activity_group_zscore.csv"), row.names = TRUE)
  write.csv(
    data.frame(pathway = names(gene_count), n_genes_used = as.integer(gene_count), row.names = NULL),
    file.path(tab_dir, "fig3H_pathway_activity_gene_counts.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      parameter = c("residualize_nCount", "androgen_fallback", "z_clip_quantile", "z_clip_abs"),
      value = c(use_ncount_resid, use_androgen_fallback, z_clip_q, z_abs),
      row.names = NULL
    ),
    file.path(tab_dir, "fig3H_pathway_activity_settings.csv"),
    row.names = FALSE
  )

  # Match paper-style panel: compact heatmap, left/top dendrograms, right row labels,
  # horizontal color key at top-left, and panel letter H.
  col_fun <- colorRamp2(c(-z_abs, 0, z_abs), c("#38807F", "#f0f0f0", "#AB2329"))

  # Keep paper-like visible order while retaining clustering branches.
  row_pref <- c(
    "NFkB", "VEGF", "TNFa", "MAPK", "EGFR", "JAK-STAT", "p53",
    "TGFb", "Androgen", "WNT", "PI3K", "Estrogen", "Hypoxia", "Trail"
  )
  row_pref <- row_pref[row_pref %in% rownames(z_draw)]
  if (length(row_pref) > 0) {
    z_draw <- z_draw[c(row_pref, setdiff(rownames(z_draw), row_pref)), , drop = FALSE]
  }

  row_hc <- hclust(dist(z_draw))
  row_dend <- as.dendrogram(row_hc)
  row_w <- setNames(seq_along(rownames(z_draw)), rownames(z_draw))
  row_dend <- reorder(row_dend, wts = row_w, agglo.FUN = mean)

  col_hc <- hclust(dist(t(z_draw)))
  col_dend <- as.dendrogram(col_hc)
  col_w <- setNames(c(1, 2, 3), c("DTStem", "LowStem", "HighStem"))
  col_w <- col_w[colnames(z_draw)]
  col_dend <- reorder(col_dend, wts = col_w, agglo.FUN = mean)

  ht <- Heatmap(
    z_draw,
    name = "Value",
    col = col_fun,
    cluster_rows = row_dend,
    cluster_columns = col_dend,
    row_names_side = "right",
    row_dend_side = "left",
    column_dend_side = "top",
    show_heatmap_legend = FALSE,
    row_names_gp = grid::gpar(fontsize = 14),
    column_names_gp = grid::gpar(fontsize = 14),
    column_names_rot = 45,
    rect_gp = grid::gpar(col = NA),
    width = grid::unit(4.6, "cm"),
    height = grid::unit(7.2, "cm")
  )

  lgd <- ComplexHeatmap::Legend(
    title = "Color Key",
    col_fun = col_fun,
    at = c(-z_abs, z_abs),
    labels = c(sprintf("%.1f", -z_abs), sprintf("%.1f", z_abs)),
    direction = "horizontal",
    legend_width = grid::unit(2.0, "cm"),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9)
  )

  draw_h_panel <- function(device = c("pdf", "png")) {
    device <- match.arg(device)
    if (device == "pdf") {
      pdf(file.path(fig_dir, "Fig3H_pathway_activity_heatmap.pdf"), width = 7.2, height = 8.2, bg = "#EBEBEB")
    } else {
      png(file.path(fig_dir, "Fig3H_pathway_activity_heatmap.png"), width = 1800, height = 2050, res = 300, bg = "#EBEBEB")
    }
    grid::grid.newpage()
    draw(
      ht,
      newpage = FALSE,
      padding = grid::unit(c(15, 20, 12, 16), "mm")
    )
    draw(lgd, x = grid::unit(0.17, "npc"), y = grid::unit(0.95, "npc"), just = c("left", "top"))
    grid::grid.text(
      "H",
      x = grid::unit(0.04, "npc"),
      y = grid::unit(0.97, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontfamily = "serif", fontface = "bold", fontsize = 40)
    )
    dev.off()
  }

  draw_h_panel("pdf")
  draw_h_panel("png")
  TRUE
  })

  # 17 Fig3I: metabolic pathway bubble (group-level activity)
  restore_or_run(17, "fig3I_metabolic_bubble", {
  expr <- mal_data$expr
  meta <- mal_data$seu@meta.data
  groups <- as.character(meta$stem_group)
  names(groups) <- rownames(meta)
  common_cells <- intersect(colnames(expr), names(groups))
  expr <- expr[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  groups <- groups[common_cells]
  use_ncount_resid_i <- tolower(Sys.getenv("FIG3I_RESIDUALIZE_NCOUNT", unset = "1")) %in% c("1", "true", "yes")
  ncount_cov_i <- NULL
  if (use_ncount_resid_i && "nCount_RNA" %in% colnames(meta)) {
    ncount_cov_i <- log1p(as.numeric(meta$nCount_RNA))
  }

  klist <- NULL
  if (!is.null(kegg_pathway_name_map) && length(kegg_pathway_name_map) > 0) {
    klist <- kegg_pathway_name_map
  } else {
    klist <- tryCatch(keggList("pathway", "hsa"), error = function(e) NULL)
  }
  if (!is.null(klist) && length(klist) > 0) {
    kdf <- data.frame(
      path_id = names(klist),
      path_name = as.character(klist),
      stringsAsFactors = FALSE
    )
    kdf$path_name <- sub(" - Homo sapiens \\(human\\)$", "", kdf$path_name)
    met_pat <- "(?i)metabolism|biosynthesis|degradation|glycolysis|gluconeogenesis|citrate cycle|fatty acid|amino acid|glutathione|steroid|pyruvate|butanoate|propanoate"
    kdf <- kdf[grepl(met_pat, kdf$path_name, perl = TRUE), , drop = FALSE]
  } else if (!is.null(local_kegg_map) && length(local_kegg_map) > 0) {
    kdf <- data.frame(
      path_id = names(local_kegg_map),
      path_name = names(local_kegg_map),
      stringsAsFactors = FALSE
    )
    # Offline-safe naming: avoid per-pathway KEGG online requests when DNS/network is unstable.
    if (!is.null(kegg_pathway_name_map) && length(kegg_pathway_name_map) > 0) {
      idx <- match(kdf$path_id, names(kegg_pathway_name_map))
      hit <- !is.na(idx)
      kdf$path_name[hit] <- unname(kegg_pathway_name_map[idx[hit]])
    }
    kdf$path_name[is.na(kdf$path_name) | !nzchar(kdf$path_name)] <- kdf$path_id
    pid_num <- suppressWarnings(as.integer(sub("^hsa", "", kdf$path_id)))
    kdf <- kdf[!is.na(pid_num) & pid_num < 2000, , drop = FALSE]
    message("Fig3I: KEGG list unavailable, using offline hsa<2000 metabolism proxy from scPagwas.")
  } else {
    stop("Fig3I: cannot load KEGG pathways (online and local both unavailable).")
  }

  score_rows <- list()
  debug_meta <- list()
  for (i in seq_len(nrow(kdf))) {
    pid <- kdf$path_id[i]
    nm <- as.character(kdf$path_name[i])
    if (!nzchar(nm) || grepl("^hsa\\d{5}$", nm)) {
      nm <- get_kegg_name(pid)
    }
    nm <- sub(" - Homo sapiens \\(human\\)$", "", nm)
    nm <- sub(" \\- reference pathway$", "", nm)
    syms <- get_kegg_symbols(pid)
    overlap <- length(intersect(toupper(rownames(expr)), syms))
    debug_meta[[length(debug_meta) + 1]] <- data.frame(
      pathway = nm, kegg_id = pid, n_symbols = length(syms), n_overlap = overlap, stringsAsFactors = FALSE
    )
    vals <- score_pathway_by_cells(expr, syms)
    if (is.null(vals)) next
    if (!is.null(ncount_cov_i)) {
      vals <- residualize_vector(as.numeric(vals), ncount_cov_i, min_n = 30L)
      names(vals) <- common_cells
    }
    tmp <- data.frame(
      cell_id = names(vals),
      pathway_id = pid,
      pathway = nm,
      score = as.numeric(vals),
      stringsAsFactors = FALSE
    )
    score_rows[[length(score_rows) + 1]] <- tmp
  }
  write.csv(bind_rows(debug_meta), file.path(tab_dir, "fig3I_kegg_overlap_debug.csv"), row.names = FALSE)
  score_long <- bind_rows(score_rows)
  if (nrow(score_long) == 0) {
    stop("Fig3I: no metabolic pathways scored. See fig3I_kegg_overlap_debug.csv")
  }
  score_long$stem_group <- groups[score_long$cell_id]
  score_long <- score_long[!is.na(score_long$stem_group), , drop = FALSE]

  agg <- score_long %>%
    group_by(pathway_id, pathway, stem_group) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

  group_order <- c("DTStem", "HighStem", "LowStem")
  agg <- agg %>%
    mutate(stem_group = factor(stem_group, levels = group_order))

  # Rank pathways by HighStem advantage.
  rank_df <- agg %>%
    tidyr::pivot_wider(names_from = stem_group, values_from = score) %>%
    mutate(high_adv = HighStem - rowMeans(cbind(DTStem, LowStem), na.rm = TRUE)) %>%
    arrange(desc(high_adv))

  top_n <- suppressWarnings(as.integer(Sys.getenv("FIG3I_TOP_N", unset = "90")))
  if (is.na(top_n) || top_n < 20) top_n <- 90L
  keep_ids <- head(rank_df$pathway_id, n = min(top_n, nrow(rank_df)))
  agg <- agg %>% filter(pathway_id %in% keep_ids)

  # Per-pathway z-score then map to [0,1] via normal CDF for richer group contrasts.
  norm_df <- agg %>%
    group_by(pathway_id, pathway) %>%
    mutate(
      mu = mean(score, na.rm = TRUE),
      sdv = sd(score, na.rm = TRUE),
      z = ifelse(is.na(sdv) | sdv < 1e-12, 0, (score - mu) / sdv),
      value = pnorm(z)
    ) %>%
    ungroup()

  # Keep ranking order in y-axis.
  path_order <- rev(rank_df$pathway[match(keep_ids, rank_df$pathway_id)])
  norm_df$pathway <- factor(norm_df$pathway, levels = path_order)

  write.csv(norm_df, file.path(tab_dir, "fig3I_metabolic_bubble_values.csv"), row.names = FALSE)

  p <- ggplot(norm_df, aes(x = stem_group, y = pathway, size = value, color = value)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(
      range = c(1.0, 7.2),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) +
    scale_color_gradientn(
      colors = c("#67C5CF", "#F0E442", "#E53935"),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) +
    labs(x = NULL, y = NULL, size = "Value", color = "Value", tag = "I") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.x = element_line(color = "#e6e6e6", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.8),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 11, angle = 65, hjust = 1, vjust = 1),
      plot.tag = element_text(family = "serif", face = "bold", size = 34),
      plot.tag.position = c(0.01, 0.99),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    )
  save_plot_both(p, "Fig3I_metabolic_pathway_bubble", 8.8, 13)
  TRUE
  })
} else {
  cat("ℹ skip_hi=ON: skipping Fig3H/Fig3I malignant-only pathway/metabolism steps.\n")
}

# 18 provenance summary
restore_or_run(18, "fig3_provenance", {
  lines <- c(
    "Fig3 analysis provenance",
    paste0("Run dir: ", run_dir),
    "Methods lock from paper: CellChat (ligand-receptor inference), netVisual_circle/netVisual_bubble, outgoing/incoming role analysis.",
    "This run follows method flow, not paper result matching.",
    paste0("CellChat labels: ", paste(levels(cellchat@idents), collapse = ", ")),
    paste0("Interactions exported: ", file.path(tab_dir, "fig3_cellchat_interactions_all.csv")),
    "Fig3H method: KEGG-defined oncogenic pathway gene sets + per-cell mean expression score + group aggregation + z-score heatmap.",
    "Fig3I method: KEGG metabolic pathways + per-cell mean expression score + group aggregation + per-pathway z-score then pnorm() mapping to [0,1].",
    "Fig3C note: outgoing/incoming strengths are computed from CellChat weight matrix after removing self-loops and clipping only top 5% tails for readability.",
    "Fig3F/G note: use significant HighStem-related LR pairs and rank by communication probability for compact bubble panels."
  )
  writeLines(lines, con = file.path(tab_dir, "fig3_analysis_provenance.txt"))
  TRUE
})

cat("✅ Fig3 method-flow analysis complete.\n")
cat("Figures dir:", fig_dir, "\n")
cat("Tables dir:", tab_dir, "\n")
