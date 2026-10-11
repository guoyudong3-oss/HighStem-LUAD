# Fig S3 A, B: redrawn with CellChat functions from the main CellChat object (checkpoint 61); analysis parameters unchanged
suppressPackageStartupMessages({library(CellChat); library(ComplexHeatmap); library(grid)})
V <- "results/reproduce_paper/runs/final_malignant_v3_20261009/"
source("03_cellchat/_get_colors_extract.R")
cellchat <- readRDS(paste0(V, "J4/main_1200/checkpoints_fig3/61_cellchat_compute_centrality.rds"))
out <- paste0(V, "W/FigS3/panels/")
lev <- levels(cellchat@idents); cols <- get_fig3_group_colors(lev); gs <- as.numeric(table(cellchat@idents))
short <- c("B lymphocytes" = "B", "DTStem" = "DTStem", "Endothelial cells" = "Endo", "Fibroblasts" = "Fib", "HighStem" = "HighStem", "LowStem" = "LowStem", "MAST cells" = "Mast",
           "Myeloid cells" = "Myeloid", "NK cells" = "NK", "NormalEpi" = "NormalEpi", "Oligodendrocytes" = "Oligodend.", "T lymphocytes" = "T", "T/NK cells" = "T/NK", "Undetermined" = "Undeterm.")
stopifnot(all(lev %in% names(short)))
n <- length(lev); ang <- 2 * pi * (seq_len(n) - 1) / n; cx <- cos(ang); cy <- sin(ang)
nx <- (cx - min(cx)) / (max(cx) - min(cx)) * 2 - 1; ny <- (cy - min(cy)) / (max(cy) - min(cy)) * 2 - 1
# copy of netVisual_circle with igraph vertex labels switched off; node names are added outside the circle
nvc_src <- deparse(netVisual_circle)
k <- grep("^    plot\\(g, ", nvc_src); stopifnot(length(k) == 1)
nvc_src[k] <- sub("plot\\(g, ", "plot(g, vertex.label = NA, ", nvc_src[k])
netVisual_circle_nolabel <- eval(parse(text = nvc_src)); environment(netVisual_circle_nolabel) <- environment(netVisual_circle)
draw_circle <- function(mat, title) {
  netVisual_circle_nolabel(mat, color.use = cols, vertex.weight = gs, weight.scale = TRUE, label.edge = FALSE, title.name = NULL, margin = 0.34)
  for (i in seq_len(n)) {
    r <- 1.22 + 0.14 * (1 - abs(sin(ang[i]))); adjx <- if (cx[i] < -0.05) 1 else if (cx[i] > 0.05) 0 else 0.5
    adjy <- if (sin(ang[i]) > 0.9) 0 else if (sin(ang[i]) < -0.9) 1 else 0.5
    if (abs(sin(ang[i])) < 0.25) { text(nx[i] + sign(cx[i]) * 0.34, ny[i] + 0.16, short[lev[i]], adj = c(0.5, 0), cex = 1, xpd = TRUE); next }
    text(nx[i] * r, ny[i] * r, short[lev[i]], adj = c(adjx, adjy), cex = 1, xpd = TRUE)
  }
  mtext(title, side = 3, line = 0.1, cex = 1, font = 2)
}
cairo_pdf(paste0(out, "S3A_circles_v2.pdf"), width = 104 / 25.4, height = 60 / 25.4, pointsize = 8, family = "Arial")
par(mfrow = c(1, 2), xpd = TRUE, mar = c(0.3, 2.0, 1.4, 2.4))
draw_circle(cellchat@net$count, "Number of interactions")
draw_circle(cellchat@net$weight, "Interaction strength")
dev.off()
# B: signalling-role heatmaps; top 20 pathways ranked by summed netP probability (same ranking as 18b)
ps <- sort(apply(cellchat@netP$prob, 3, sum, na.rm = TRUE), decreasing = TRUE)
sig <- names(ps)[seq_len(20)]
write.csv(data.frame(rank = seq_along(ps), pathway = names(ps), prob_sum = as.numeric(ps)), paste0(out, "S3B_pathway_rank_all_v2.csv"), row.names = FALSE)
# copy of netAnalysis_signalingRole_heatmap that only shortens group names
hm_src <- deparse(netAnalysis_signalingRole_heatmap)
k <- grep("dimnames\\(outgoing\\) <- list\\(levels\\(object@idents\\), names\\(centr\\)\\)", hm_src); stopifnot(length(k) == 1)
hm_src[k] <- "    dimnames(outgoing) <- list(unname(short[levels(object@idents)]), names(centr))"
k2 <- grep("ha2 = HeatmapAnnotation\\(Strength = anno_barplot\\(colSums\\(mat.ori\\), ", hm_src); stopifnot(length(k2) == 1)
hm_src[k2] <- sub("anno_barplot\\(colSums\\(mat.ori\\), ", "anno_barplot(colSums(mat.ori), height = grid::unit(7, 'mm'), ", hm_src[k2])
role_heatmap_short <- eval(parse(text = hm_src)); environment(role_heatmap_short) <- environment(netAnalysis_signalingRole_heatmap)
cols2 <- setNames(cols, unname(short[names(cols)]))
mk <- function(p) role_heatmap_short(cellchat, signaling = sig, pattern = p, color.use = unname(cols2), width = 4.2, height = 5.2, font.size = 8, font.size.title = 8)
ho <- mk("outgoing"); hi <- mk("incoming")
cairo_pdf(paste0(out, "S3B_role_heatmaps_v2.pdf"), width = 180 / 25.4, height = 88 / 25.4, pointsize = 8, family = "Arial")
draw(ho + hi, ht_gap = unit(14, "mm"), padding = unit(c(2, 6, 4, 2), "mm"))
dev.off()
# export totals and pathway ranking for cross-checking
w <- cellchat@net$weight; cnt <- cellchat@net$count
write.csv(data.frame(group = lev, n_cells = gs, out_count = rowSums(cnt), in_count = colSums(cnt), out_weight = rowSums(w), in_weight = colSums(w)), paste0(out, "S3A_group_totals_v2.csv"), row.names = FALSE)
cat("pathways:", sig, "\n")
