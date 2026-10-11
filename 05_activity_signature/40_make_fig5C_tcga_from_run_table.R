#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggpubr)
})

get_arg <- function(args, key, default = NA_character_) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", key, "="), "", hit[[1]]))
  idx <- which(args == key)
  if (length(idx) > 0 && idx[[1]] < length(args)) return(args[idx[[1]] + 1])
  default
}

args <- commandArgs(trailingOnly = TRUE)
score_csv <- get_arg(args, "--score-csv", "")
out_png <- get_arg(args, "--out-png", "")
out_pdf <- get_arg(args, "--out-pdf", "")
out_summary <- get_arg(args, "--out-summary", "")
title_txt <- get_arg(args, "--title", "TCGA-LUAD: HighStem activity")

if (!nzchar(score_csv) || !file.exists(score_csv)) stop("Missing --score-csv")
if (!nzchar(out_png)) stop("Missing --out-png")
if (!nzchar(out_pdf)) stop("Missing --out-pdf")
if (!nzchar(out_summary)) stop("Missing --out-summary")

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)

dt <- fread(score_csv)
if (!all(c("PlasticStem_Score", "sample_type") %in% names(dt))) {
  stop("score csv missing PlasticStem_Score / sample_type")
}

dt <- dt[is.finite(PlasticStem_Score)]
dt[, group := fifelse(tolower(sample_type) %in% c("normal", "solid tissue normal"), "Normal", "LUAD")]
dt <- dt[group %in% c("Normal", "LUAD")]
dt[, group := factor(group, levels = c("Normal", "LUAD"))]

pw <- wilcox.test(PlasticStem_Score ~ group, data = dt)$p.value

p <- ggplot(dt, aes(x = group, y = PlasticStem_Score, fill = group, color = group)) +
  geom_violin(width = 0.85, alpha = 0.35, linewidth = 0) +
  geom_boxplot(width = 0.10, outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_jitter(width = 0.08, size = 0.45, alpha = 0.85, show.legend = FALSE) +
  stat_compare_means(comparisons = list(c("Normal", "LUAD")), label = "p.signif", method = "wilcox.test", size = 4) +
  scale_fill_manual(values = c("Normal" = "#7CC7C9", "LUAD" = "#D96B86")) +
  scale_color_manual(values = c("Normal" = "#0C9BA4", "LUAD" = "#DC1F3A")) +
  labs(title = title_txt, x = NULL, y = "HighStem") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

ggsave(out_pdf, p, width = 202 / 72, height = 282 / 72, units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
ggsave(out_png, p, width = 202 / 72, height = 282 / 72, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

fwrite(
  data.table(
    n_total = nrow(dt),
    n_normal = sum(dt$group == "Normal"),
    n_luad = sum(dt$group == "LUAD"),
    wilcox_p = pw
  ),
  out_summary
)

cat("OUT_PNG=", normalizePath(out_png, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("OUT_PDF=", normalizePath(out_pdf, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("OUT_SUMMARY=", normalizePath(out_summary, winslash = "/", mustWork = TRUE), "\n", sep = "")
