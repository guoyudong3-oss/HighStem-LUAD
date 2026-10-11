#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(infercnv)
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

fit_gmm2 <- function(x, max_iter = 200, tol = 1e-7) {
  x <- x[is.finite(x)]
  if (length(x) < 100) stop("Too few values for GMM: ", length(x))
  q <- quantile(x, c(0.3, 0.7), na.rm = TRUE)
  mu <- as.numeric(q)
  sdv <- rep(stats::sd(x), 2)
  pi <- c(0.5, 0.5)
  ll_prev <- -Inf
  for (iter in seq_len(max_iter)) {
    d1 <- pi[1] * dnorm(x, mean = mu[1], sd = sdv[1])
    d2 <- pi[2] * dnorm(x, mean = mu[2], sd = sdv[2])
    den <- d1 + d2 + 1e-300
    r1 <- d1 / den
    r2 <- 1 - r1
    n1 <- sum(r1); n2 <- sum(r2)
    pi <- c(n1, n2) / length(x)
    mu <- c(sum(r1 * x) / n1, sum(r2 * x) / n2)
    v1 <- sum(r1 * (x - mu[1])^2) / n1
    v2 <- sum(r2 * (x - mu[2])^2) / n2
    sdv <- sqrt(pmax(c(v1, v2), 1e-12))
    ll <- sum(log(den))
    if (abs(ll - ll_prev) < tol) break
    ll_prev <- ll
  }
  list(pi = pi, mu = mu, sd = sdv)
}

gmm_posterior_cutoff <- function(x) {
  gmm <- fit_gmm2(x)
  o <- order(gmm$mu)
  pi1 <- gmm$pi[o[1]]; mu1 <- gmm$mu[o[1]]; sd1 <- gmm$sd[o[1]]
  pi2 <- gmm$pi[o[2]]; mu2 <- gmm$mu[o[2]]; sd2 <- gmm$sd[o[2]]
  f <- function(z) (log(pi1) - log(sd1) - (z - mu1)^2 / (2 * sd1^2)) - (log(pi2) - log(sd2) - (z - mu2)^2 / (2 * sd2^2))
  lo <- min(mu1, mu2); hi <- max(mu1, mu2)
  if (is.finite(f(lo)) && is.finite(f(hi)) && f(lo) * f(hi) < 0) {
    return(uniroot(f, lower = lo, upper = hi)$root)
  }
  NA_real_
}

args <- commandArgs(trailingOnly = TRUE)
run_dir <- get_arg(args, "--run-dir", NULL)
method <- get_arg(args, "--method", "ref_mean_3sd")

if (is.null(run_dir) || !nzchar(run_dir)) {
  stop("Usage: Rscript 00c5_export_malignant_set_from_runfinal.R --run-dir <path> [--method <name>]\n",
       "Methods: all_epi_median_1.5mad, all_epi_median_3mad, ref_q99, ref_mean_3sd, obs_gmm2, obs_top_5pct, obs_top_10pct")
}

run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
infer_dir <- file.path(run_dir, "infercnv_output")
tables_dir <- file.path(run_dir, "tables")
rds_dir <- file.path(run_dir, "rds")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

ann_file <- file.path(tables_dir, "infercnv_annotations_paper_exact.tsv")
meta0_file <- file.path(tables_dir, "cnv_paper_exact_cell_metadata.csv")
infer_rds <- file.path(infer_dir, "run.final.infercnv_obj")
seu_in <- file.path(rds_dir, "epithelial_plus_refs_for_infercnv_paper_exact.rds")

if (!file.exists(ann_file)) stop("Missing: ", ann_file)
if (!file.exists(meta0_file)) stop("Missing: ", meta0_file)
if (!file.exists(infer_rds)) stop("Missing: ", infer_rds)

ann <- read.table(ann_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(ann) <- c("cell_id", "infercnv_group")

meta0 <- read.csv(meta0_file, stringsAsFactors = FALSE)
cell_type_map <- meta0$Cell_type; names(cell_type_map) <- meta0$cell_id
origin_map <- meta0$Sample_Origin; names(origin_map) <- meta0$cell_id

infer_obj <- readRDS(infer_rds)
expr <- infer_obj@expr.data
if (!is.matrix(expr)) stop("infer_obj@expr.data is not a matrix.")
cnv_score <- colMeans(abs(expr - 1))

ann$cnv_score_runfinal <- cnv_score[match(ann$cell_id, names(cnv_score))]
ann$Cell_type <- unname(cell_type_map[ann$cell_id])
ann$Sample_Origin <- unname(origin_map[ann$cell_id])

ref_scores <- ann$cnv_score_runfinal[ann$infercnv_group == "Ref_Epithelial"]
obs_scores <- ann$cnv_score_runfinal[ann$infercnv_group == "Observation"]
if (length(ref_scores) < 10) stop("Too few reference epithelial cells: ", length(ref_scores))
if (length(obs_scores) < 10) stop("Too few observation epithelial cells: ", length(obs_scores))

threshold <- NA_real_
threshold_basis <- ""
method <- tolower(method)

if (method == "all_epi_median_1.5mad") {
  threshold <- median(ann$cnv_score_runfinal, na.rm = TRUE) + 1.5 * mad(ann$cnv_score_runfinal, na.rm = TRUE)
  threshold_basis <- "all epithelial (ref+obs): median + 1.5*MAD"
} else if (method == "all_epi_median_3mad") {
  threshold <- median(ann$cnv_score_runfinal, na.rm = TRUE) + 3.0 * mad(ann$cnv_score_runfinal, na.rm = TRUE)
  threshold_basis <- "all epithelial (ref+obs): median + 3*MAD"
} else if (method == "ref_q99") {
  threshold <- as.numeric(quantile(ref_scores, 0.99, na.rm = TRUE))
  threshold_basis <- "reference epithelial: 99th percentile"
} else if (method == "ref_mean_3sd") {
  threshold <- mean(ref_scores, na.rm = TRUE) + 3.0 * sd(ref_scores, na.rm = TRUE)
  threshold_basis <- "reference epithelial: mean + 3*SD"
} else if (method == "obs_gmm2") {
  threshold <- gmm_posterior_cutoff(obs_scores)
  threshold_basis <- "observation epithelial: 2-Gaussian mixture posterior=0.5 cutoff"
} else if (method == "obs_top_5pct") {
  threshold <- as.numeric(quantile(obs_scores, 0.95, na.rm = TRUE))
  threshold_basis <- "observation epithelial: top 5% CNV score"
} else if (method == "obs_top_10pct") {
  threshold <- as.numeric(quantile(obs_scores, 0.90, na.rm = TRUE))
  threshold_basis <- "observation epithelial: top 10% CNV score"
} else {
  stop("Unknown --method: ", method)
}

if (!is.finite(threshold)) stop("Computed threshold is not finite for method: ", method)

ann$malignant_status <- ifelse(
  ann$infercnv_group == "Observation" & ann$cnv_score_runfinal > threshold,
  "malignant",
  "normal_epithelial"
)

out_meta <- file.path(tables_dir, paste0("cnv_paper_exact_cell_metadata_runfinal_", method, ".csv"))
out_obs <- file.path(tables_dir, paste0("cnv_paper_exact_observation_scores_runfinal_", method, ".csv"))
out_summary <- file.path(tables_dir, paste0("cnv_paper_exact_summary_runfinal_", method, ".txt"))
out_seu <- file.path(rds_dir, paste0("malignant_epithelial_paper_exact_runfinal_", method, ".rds"))

meta_out <- ann[, c("cell_id", "Cell_type", "Sample_Origin", "infercnv_group", "cnv_score_runfinal", "malignant_status")]
write.csv(meta_out, out_meta, row.names = FALSE)

obs_out <- meta_out[meta_out$infercnv_group == "Observation", c("cell_id", "Sample_Origin", "cnv_score_runfinal", "malignant_status")]
obs_out <- obs_out[order(-obs_out$cnv_score_runfinal), ]
write.csv(obs_out, out_obs, row.names = FALSE)

seu_msg <- "Skipped Seurat export."
if (file.exists(seu_in)) {
  suppressPackageStartupMessages(library(Seurat))
  seu <- readRDS(seu_in)
  mal_cells <- meta_out$cell_id[meta_out$malignant_status == "malignant"]
  mal_cells <- intersect(mal_cells, colnames(seu))
  mal_obj <- subset(seu, cells = mal_cells)
  saveRDS(mal_obj, out_seu)
  seu_msg <- paste0("Saved malignant Seurat RDS: ", out_seu, " (cells=", length(mal_cells), ")")
}

capture.output({
  cat("inferCNV malignant export from run.final.infercnv_obj\n")
  cat("run_dir:", run_dir, "\n")
  cat("method:", method, "\n")
  cat("threshold basis:", threshold_basis, "\n")
  cat("threshold:", threshold, "\n\n")
  cat("status counts:\n")
  print(table(meta_out$malignant_status, useNA = "ifany"))
  cat("\ninfercnv groups:\n")
  print(table(meta_out$infercnv_group, useNA = "ifany"))
  cat("\nObservation malignant:\n")
  cat(sum(meta_out$infercnv_group == "Observation" & meta_out$malignant_status == "malignant", na.rm = TRUE), "\n")
  cat("\nscore summary (all epithelial):\n")
  print(summary(meta_out$cnv_score_runfinal))
  cat("\nreference score summary:\n")
  print(summary(ref_scores))
  cat("\nobservation score summary:\n")
  print(summary(obs_scores))
  cat("\n", seu_msg, "\n", sep = "")
}, file = out_summary)

cat("Saved:", out_meta, "\n")
cat("Saved:", out_obs, "\n")
cat("Saved:", out_summary, "\n")
cat(seu_msg, "\n")

