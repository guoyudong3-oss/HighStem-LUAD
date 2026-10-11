# Fig3D: pathway activity heatmap (fig3H_pathway_activity_group_zscore.csv; same clustering as 18b)
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize); library(grid)})
V <- "results/reproduce_paper/runs/final_malignant_v3_20261009/"
z <- as.matrix(read.csv(paste0(V, "J4/main_1200/tables/fig3H_pathway_activity_group_zscore.csv"), row.names = 1, check.names = FALSE))
z_draw <- z[, c("DTStem", "LowStem", "HighStem")]; z_abs <- max(abs(z_draw)); z_abs <- ceiling(z_abs * 10) / 10
col_fun <- colorRamp2(c(-z_abs, 0, z_abs), c("#38807F", "#f0f0f0", "#AB2329"))
row_pref <- c("NFkB", "VEGF", "TNFa", "MAPK", "EGFR", "JAK-STAT", "p53", "TGFb", "Androgen", "WNT", "PI3K", "Estrogen", "Hypoxia", "Trail")
z_draw <- z_draw[c(row_pref[row_pref %in% rownames(z_draw)], setdiff(rownames(z_draw), row_pref)), , drop = FALSE]
row_dend <- reorder(as.dendrogram(hclust(dist(z_draw))), wts = setNames(seq_len(nrow(z_draw)), rownames(z_draw)), agglo.FUN = mean)
col_w <- setNames(c(1, 2, 3), c("DTStem", "LowStem", "HighStem"))[colnames(z_draw)]
col_dend <- reorder(as.dendrogram(hclust(dist(t(z_draw)))), wts = col_w, agglo.FUN = mean)
ht <- Heatmap(z_draw, name = "z", col = col_fun, cluster_rows = row_dend, cluster_columns = col_dend, row_names_side = "right", row_dend_side = "left", column_dend_side = "top",
  show_heatmap_legend = FALSE, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 0.5),
  width = unit(24, "mm"), height = unit(46, "mm"), row_dend_width = unit(10, "mm"), column_dend_height = unit(7, "mm"))
lgd <- Legend(title = "Standardized score", col_fun = col_fun, at = c(-z_abs, 0, z_abs), labels = c(sprintf("%.1f", -z_abs), "0", sprintf("%.1f", z_abs)), direction = "horizontal",
  legend_width = unit(26, "mm"), title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8))
out <- paste0(V, "W/Fig3/")
cairo_pdf(paste0(out, "Fig3D_heatmap.pdf"), width = 80 / 25.4, height = 88 / 25.4, family = "Arial")
grid.newpage(); draw(ht, newpage = FALSE, padding = unit(c(24, 2, 4, 2), "mm")); draw(lgd, x = unit(14, "mm"), y = unit(12.5, "mm"), just = c("left", "top"))
dev.off()
