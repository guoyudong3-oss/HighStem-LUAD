#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
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

args <- commandArgs(trailingOnly = TRUE)

core_genes <- c("CST6", "MUC5B", "TIMP1", "S100A16", "MARCKSL1")
axis_genes <- c("PDPK1", "AKT1", "MTOR")
summary_genes <- c("PIK3CA", "PIK3CB", "PDPK1", "AKT1", "MTOR")
focus_pathways <- c("PI3K", "Hypoxia", "NFkB", "JAK-STAT")

expr_file <- get_arg(args, "--expr-file", "data/external/tcga_luad/expression.tsv.gz")
mapping_file <- get_arg(args, "--mapping-file", "data/external/ensembl_to_symbol.csv")
pathway_z_file <- get_arg(args, "--pathway-z-file", "results/reproduce_paper/runs/paper_main_tme_cellchat_highratio10_rawlog1p_allraw_20260215_013749/tables/fig3H_pathway_activity_group_zscore.csv")
pathway_mean_file <- get_arg(args, "--pathway-mean-file", "results/reproduce_paper/runs/paper_main_tme_cellchat_highratio10_rawlog1p_allraw_20260215_013749/tables/fig3H_pathway_activity_group_mean.csv")
gene_file <- get_arg(args, "--gene-file", "results/reproduce_paper/tables/final_consensus_genes.csv")

fig_dir <- get_arg(args, "--fig-dir", "results/reproduce_paper/figures")
cache_dir <- get_arg(args, "--cache-dir", "results/reproduce_paper/fig8")
paper_dir <- get_arg(args, "--paper-dir", "paper_figures")
out_prefix <- get_arg(args, "--out-prefix", "Fig8_bioinfo_PI3K_allraw_full")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
if (!is.null(paper_dir) && nzchar(paper_dir)) {
  dir.create(paper_dir, recursive = TRUE, showWarnings = FALSE)
}

out_png <- file.path(fig_dir, paste0(out_prefix, ".png"))
out_pdf <- file.path(fig_dir, paste0(out_prefix, ".pdf"))
paper_png <- if (!is.null(paper_dir) && nzchar(paper_dir)) file.path(paper_dir, paste0(out_prefix, ".png")) else ""
paper_pdf <- if (!is.null(paper_dir) && nzchar(paper_dir)) file.path(paper_dir, paste0(out_prefix, ".pdf")) else ""
paper_png_canonical <- if (!is.null(paper_dir) && nzchar(paper_dir)) file.path(paper_dir, "Fig8.png") else ""
paper_pdf_canonical <- if (!is.null(paper_dir) && nzchar(paper_dir)) file.path(paper_dir, "Fig8.pdf") else ""

expr_cache_file <- file.path(cache_dir, "fig8_bioinfo_tcga_tumor_log2tpm_core_pi3k.csv")
scatter_cor_file <- file.path(cache_dir, "fig8_scatter_core_vs_pi3k_axis_spearman.csv")
score_cor_file <- file.path(cache_dir, "fig8_signature_vs_pi3k_axis_spearman.csv")
core_axis_heatmap_file <- file.path(cache_dir, "fig8_core_vs_pi3k_axis_heatmap_spearman.csv")

sample_type_from_barcode <- function(barcode) {
  x <- tryCatch(as.integer(substr(strsplit(barcode, "-")[[1]][4], 1, 2)), error = function(e) NA_integer_)
  if (is.na(x)) return("Unknown")
  if (x < 10) return("Tumor")
  if (x < 20) return("Normal")
  "Other"
}

zscore_vec <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-99) return("< 1e-99")
  formatC(p, format = "e", digits = 2)
}

load_tumor_expr <- function() {
  expr <- read.delim(gzfile(expr_file), check.names = FALSE, stringsAsFactors = FALSE)
  mapdf <- read.csv(mapping_file, stringsAsFactors = FALSE)
  mapdf$gene_id <- sub("\\..*$", "", mapdf$gene_id)

  colnames(expr)[1] <- "gene_id"
  expr$gene_id <- sub("\\..*$", "", expr$gene_id)
  idx <- match(expr$gene_id, mapdf$gene_id)
  expr$symbol <- mapdf$gene_name[idx]
  expr <- expr[!is.na(expr$symbol) & expr$symbol != "", , drop = FALSE]

  num <- as.matrix(expr[, setdiff(colnames(expr), c("gene_id", "symbol")), drop = FALSE])
  storage.mode(num) <- "numeric"

  agg <- rowsum(num, group = expr$symbol)
  counts <- table(expr$symbol)[rownames(agg)]
  agg <- agg / as.vector(counts)
  agg <- log2(agg + 1)

  tumor_keep <- vapply(colnames(agg), sample_type_from_barcode, character(1)) == "Tumor"
  agg[, tumor_keep, drop = FALSE]
}

make_scatter_plot <- function(df, x_gene, y_gene) {
  ct <- suppressWarnings(stats::cor.test(df[[x_gene]], df[[y_gene]], method = "spearman", exact = FALSE))
  label <- sprintf("rho = %.3f\np = %s", unname(ct$estimate), fmt_p(ct$p.value))

  ggplot(df, aes(x = .data[[x_gene]], y = .data[[y_gene]])) +
    geom_point(color = "#444444", alpha = 0.55, size = 1.1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, color = "#c0392b") +
    annotate("text", x = Inf, y = Inf, label = label, hjust = 1.05, vjust = 1.15,
             size = 2.55, color = "#c0392b") +
    labs(
      x = sprintf("%s Expression (log2(TPM+1))", x_gene),
      y = sprintf("%s Expression (log2(TPM+1))", y_gene)
    ) +
    theme_bw(base_size = 8.5) +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 4, 6),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5)
    )
}

make_text_title <- function(title, subtitle = NULL, big = FALSE) {
  ggplot() +
    theme_void() +
    annotate("text", x = 0, y = 1, label = title, hjust = 0, vjust = 1,
             fontface = "bold", size = if (big) 6 else 5) +
    {
      if (!is.null(subtitle)) {
        annotate("text", x = 0, y = 0.55, label = subtitle, hjust = 0, vjust = 1,
                 size = if (big) 4 else 3.6, color = "#444444")
      }
    } +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")
}

load_pathway_focus <- function() {
  zmat <- read.csv(pathway_z_file, row.names = 1, check.names = FALSE)
  mmat <- read.csv(pathway_mean_file, row.names = 1, check.names = FALSE)
  zmat <- zmat[focus_pathways, c("LowStem", "DTStem", "HighStem"), drop = FALSE]
  mmat <- mmat[focus_pathways, c("LowStem", "DTStem", "HighStem"), drop = FALSE]

  out <- data.frame(
    pathway = rep(rownames(zmat), times = ncol(zmat)),
    group = rep(colnames(zmat), each = nrow(zmat)),
    zscore = as.vector(as.matrix(zmat)),
    mean_value = as.vector(as.matrix(mmat)),
    stringsAsFactors = FALSE
  )
  out$pathway <- factor(out$pathway, levels = rev(focus_pathways))
  out$group <- factor(out$group, levels = c("LowStem", "DTStem", "HighStem"))
  out
}

make_heatmap_plot <- function(pathway_df) {
  ggplot(pathway_df, aes(x = group, y = pathway, fill = zscore)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", zscore)), size = 3.0) +
    scale_fill_gradient2(
      low = "#2b8cbe", mid = "#f7f7f7", high = "#d7301f",
      midpoint = 0, limits = c(-1.2, 1.2), name = "z-score"
    ) +
    labs(
      title = "A  HighStem-focused pathway activity",
      subtitle = "Discovery single-cell transcriptomic analysis"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "right",
      plot.margin = margin(6, 8, 4, 8)
    )
}

make_core_axis_heatmap <- function(cor_df) {
  cor_df$core_gene <- factor(cor_df$core_gene, levels = rev(core_genes))
  cor_df$axis_gene <- factor(cor_df$axis_gene, levels = summary_genes)

  ggplot(cor_df, aes(x = axis_gene, y = core_gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3.0) +
    scale_fill_gradient(low = "#f7f7f7", high = "#c0392b", limits = c(0, 0.55), name = "rho") +
    labs(
      title = "C  Core genes versus PI3K-axis genes",
      subtitle = "Tumor-only TCGA-LUAD (n = 519)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "right",
      plot.margin = margin(6, 8, 4, 8)
    )
}

make_score_plot <- function(score_df) {
  score_df <- score_df[order(score_df$rho), , drop = FALSE]
  score_df$gene <- factor(score_df$gene, levels = score_df$gene)
  score_df$label <- sprintf("rho = %.2f\np = %s", score_df$rho, vapply(score_df$p_value, fmt_p, character(1)))

  ggplot(score_df, aes(x = rho, y = gene)) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "#999999", linetype = "dashed") +
    geom_segment(aes(x = 0, xend = rho, y = gene, yend = gene), linewidth = 0.8, color = "#777777") +
    geom_point(size = 3.2, color = "#c0392b") +
    geom_text(aes(label = label), hjust = -0.05, size = 2.75, lineheight = 0.9, color = "#333333") +
    coord_cartesian(xlim = c(0, 0.68), clip = "off") +
    labs(
      title = "D  Five-gene program versus PI3K-axis genes",
      subtitle = "Tumor-only TCGA-LUAD (n = 519)",
      x = "Spearman rho",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(6, 36, 4, 8)
    )
}

cat("[INFO] Reading gene list\n")
gene_tbl <- read.csv(gene_file, stringsAsFactors = FALSE)
stopifnot(all(core_genes %in% gene_tbl$gene))

cat("[INFO] Loading tumor-only TCGA expression and extracting PI3K-axis genes\n")
tumor_expr <- load_tumor_expr()
needed_genes <- unique(c(core_genes, axis_genes, summary_genes))
missing_genes <- setdiff(needed_genes, rownames(tumor_expr))
if (length(missing_genes) > 0) {
  stop("Missing genes in tumor expression matrix: ", paste(missing_genes, collapse = ", "))
}
tumor_expr <- tumor_expr[needed_genes, , drop = FALSE]
write.csv(data.frame(sample = colnames(tumor_expr), t(tumor_expr), row.names = NULL),
          expr_cache_file, row.names = FALSE)

cat("[INFO] Building scatter plots\n")
scatter_df <- data.frame(sample = colnames(tumor_expr), t(tumor_expr), row.names = NULL, check.names = FALSE)
scatter_cor_rows <- list()
scatter_plots <- list()
for (yg in core_genes) {
  for (xg in axis_genes) {
    ct <- suppressWarnings(stats::cor.test(scatter_df[[xg]], scatter_df[[yg]], method = "spearman", exact = FALSE))
    scatter_cor_rows[[length(scatter_cor_rows) + 1]] <- data.frame(
      y_gene = yg,
      x_gene = xg,
      rho = unname(ct$estimate),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
    scatter_plots[[length(scatter_plots) + 1]] <- make_scatter_plot(scatter_df, xg, yg)
  }
}
scatter_cor_df <- do.call(rbind, scatter_cor_rows)
write.csv(scatter_cor_df, scatter_cor_file, row.names = FALSE)

scatter_grid <- wrap_plots(scatter_plots, ncol = length(axis_genes))
scatter_block <- make_text_title(
  "B  Core-gene correlations with representative PI3K-axis genes",
  "Rows: CST6/MUC5B/TIMP1/S100A16/MARCKSL1; Columns: PDPK1/AKT1/MTOR"
) / scatter_grid + plot_layout(heights = c(0.08, 1))

cat("[INFO] Building focused pathway support panels\n")
pathway_df <- load_pathway_focus()
heatmap_plot <- make_heatmap_plot(pathway_df)

sig_score <- colMeans(t(apply(tumor_expr[core_genes, , drop = FALSE], 1, zscore_vec)))
score_cor_rows <- list()
for (g in summary_genes) {
  ct <- suppressWarnings(stats::cor.test(sig_score, as.numeric(tumor_expr[g, ]), method = "spearman", exact = FALSE))
  score_cor_rows[[length(score_cor_rows) + 1]] <- data.frame(
    gene = g,
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}
score_cor_df <- do.call(rbind, score_cor_rows)
write.csv(score_cor_df, score_cor_file, row.names = FALSE)
score_plot <- make_score_plot(score_cor_df)

heatmap_cor_df <- rename_cols <- scatter_cor_df[, c("y_gene", "x_gene", "rho", "p_value")]
colnames(heatmap_cor_df) <- c("core_gene", "axis_gene", "rho", "p_value")
heatmap_cor_df$axis_gene <- factor(heatmap_cor_df$axis_gene, levels = axis_genes)

for (extra_gene in setdiff(summary_genes, axis_genes)) {
  extra_rows <- lapply(core_genes, function(cg) {
    ct <- suppressWarnings(stats::cor.test(scatter_df[[extra_gene]], scatter_df[[cg]], method = "spearman", exact = FALSE))
    data.frame(
      core_gene = cg,
      axis_gene = extra_gene,
      rho = unname(ct$estimate),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
  heatmap_cor_df <- rbind(heatmap_cor_df, do.call(rbind, extra_rows))
}
heatmap_cor_df$axis_gene <- factor(heatmap_cor_df$axis_gene, levels = summary_genes)
write.csv(heatmap_cor_df, core_axis_heatmap_file, row.names = FALSE)
core_axis_heatmap_plot <- make_core_axis_heatmap(heatmap_cor_df)

right_stack <- heatmap_plot / core_axis_heatmap_plot / score_plot + plot_layout(heights = c(0.9, 1.0, 1.0))

combined <- right_stack | scatter_block +
  plot_layout(widths = c(1.05, 2.45))

cat("[INFO] Saving figure outputs\n")
ggsave(out_png, combined, width = 18, height = 13.2, dpi = 320, bg = "white")
ggsave(out_pdf, combined, width = 18, height = 13.2, bg = "white")
if (nzchar(paper_png)) ggsave(paper_png, combined, width = 18, height = 13.2, dpi = 320, bg = "white")
if (nzchar(paper_pdf)) ggsave(paper_pdf, combined, width = 18, height = 13.2, bg = "white")
if (nzchar(paper_png_canonical)) ggsave(paper_png_canonical, combined, width = 18, height = 13.2, dpi = 320, bg = "white")
if (nzchar(paper_pdf_canonical)) ggsave(paper_pdf_canonical, combined, width = 18, height = 13.2, bg = "white")

cat("[DONE] Figure saved:\n")
cat("  ", out_png, "\n")
cat("  ", out_pdf, "\n")
if (nzchar(paper_png)) cat("  ", paper_png, "\n")
if (nzchar(paper_pdf)) cat("  ", paper_pdf, "\n")
if (nzchar(paper_png_canonical)) cat("  ", paper_png_canonical, "\n")
if (nzchar(paper_pdf_canonical)) cat("  ", paper_pdf_canonical, "\n")
