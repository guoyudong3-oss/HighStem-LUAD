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

args <- commandArgs(trailingOnly = TRUE)
run_dir <- get_arg(args, "--run-dir", NULL)
if (is.null(run_dir) || !nzchar(run_dir)) {
  stop("Usage: Rscript 00c4_malignant_calling_sweep_runfinal.R --run-dir <path>")
}

run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
infer_dir <- file.path(run_dir, "infercnv_output")
tables_dir <- file.path(run_dir, "tables")
figures_dir <- file.path(run_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ann_file <- file.path(tables_dir, "infercnv_annotations_paper_exact.tsv")
meta0_file <- file.path(tables_dir, "cnv_paper_exact_cell_metadata.csv")
infer_rds <- file.path(infer_dir, "run.final.infercnv_obj")

if (!file.exists(ann_file)) stop("Missing: ", ann_file)
if (!file.exists(meta0_file)) stop("Missing: ", meta0_file)
if (!file.exists(infer_rds)) stop("Missing: ", infer_rds)

ann <- read.table(ann_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(ann) <- c("cell_id", "infercnv_group")

meta0 <- read.csv(meta0_file, stringsAsFactors = FALSE)
origin_map <- meta0$Sample_Origin
names(origin_map) <- meta0$cell_id

infer_obj <- readRDS(infer_rds)
expr <- infer_obj@expr.data
if (!is.matrix(expr)) stop("infer_obj@expr.data is not a matrix.")

# Per-cell CNV score used in our pipeline: mean absolute deviation from baseline=1.
cnv_score <- colMeans(abs(expr - 1))

ann$cnv_score <- cnv_score[match(ann$cell_id, names(cnv_score))]
ann$Sample_Origin <- unname(origin_map[ann$cell_id])

ref_scores <- ann$cnv_score[ann$infercnv_group == "Ref_Epithelial"]
obs_scores <- ann$cnv_score[ann$infercnv_group == "Observation"]
if (length(ref_scores) < 10) stop("Too few reference epithelial cells: ", length(ref_scores))
if (length(obs_scores) < 10) stop("Too few observation epithelial cells: ", length(obs_scores))

mk_row <- function(method, thr) {
  mal <- sum(obs_scores > thr, na.rm = TRUE)
  data.frame(
    method = method,
    threshold = thr,
    obs_total = length(obs_scores),
    obs_malignant = mal,
    obs_malignant_pct = 100 * mal / length(obs_scores),
    stringsAsFactors = FALSE
  )
}

thresholds <- list(
  # Conservative, but common in internal pipelines: threshold from all epithelial (ref+obs).
  all_epi_median_1.5mad = median(ann$cnv_score, na.rm = TRUE) + 1.5 * mad(ann$cnv_score, na.rm = TRUE),
  all_epi_median_3mad   = median(ann$cnv_score, na.rm = TRUE) + 3.0 * mad(ann$cnv_score, na.rm = TRUE),

  # More common in papers/forums when a clean reference exists: threshold from reference distribution.
  ref_median_3mad = median(ref_scores, na.rm = TRUE) + 3.0 * mad(ref_scores, na.rm = TRUE),
  ref_q95  = as.numeric(quantile(ref_scores, 0.95, na.rm = TRUE)),
  ref_q975 = as.numeric(quantile(ref_scores, 0.975, na.rm = TRUE)),
  ref_q99  = as.numeric(quantile(ref_scores, 0.99, na.rm = TRUE)),
  ref_q995 = as.numeric(quantile(ref_scores, 0.995, na.rm = TRUE)),
  ref_mean_2sd = mean(ref_scores, na.rm = TRUE) + 2 * sd(ref_scores, na.rm = TRUE),
  ref_mean_3sd = mean(ref_scores, na.rm = TRUE) + 3 * sd(ref_scores, na.rm = TRUE)
)

# Rank-based cutoffs sometimes used in practice when reference is uncertain (sensitivity analysis).
thresholds$obs_top_1pct  <- as.numeric(quantile(obs_scores, 0.99, na.rm = TRUE))
thresholds$obs_top_5pct  <- as.numeric(quantile(obs_scores, 0.95, na.rm = TRUE))
thresholds$obs_top_10pct <- as.numeric(quantile(obs_scores, 0.90, na.rm = TRUE))

# Two-component Gaussian mixture on observation scores (a common statistical alternative).
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

gmm_cutoff <- NA_real_
gmm <- tryCatch(fit_gmm2(obs_scores), error = function(e) NULL)
if (!is.null(gmm)) {
  # Find the point where posterior is 0.5 (weighted densities equal) between the two means.
  o <- order(gmm$mu)
  pi1 <- gmm$pi[o[1]]; mu1 <- gmm$mu[o[1]]; sd1 <- gmm$sd[o[1]]
  pi2 <- gmm$pi[o[2]]; mu2 <- gmm$mu[o[2]]; sd2 <- gmm$sd[o[2]]
  f <- function(z) (log(pi1) - log(sd1) - (z - mu1)^2 / (2 * sd1^2)) - (log(pi2) - log(sd2) - (z - mu2)^2 / (2 * sd2^2))
  # uniroot requires a sign change; if not, leave NA.
  lo <- min(mu1, mu2); hi <- max(mu1, mu2)
  if (is.finite(f(lo)) && is.finite(f(hi)) && f(lo) * f(hi) < 0) {
    gmm_cutoff <- tryCatch(uniroot(f, lower = lo, upper = hi)$root, error = function(e) NA_real_)
  }
}
if (is.finite(gmm_cutoff)) thresholds$obs_gmm2_cutoff <- gmm_cutoff

sweep <- do.call(rbind, lapply(names(thresholds), function(nm) mk_row(nm, thresholds[[nm]])))
sweep <- sweep[order(sweep$obs_malignant), ]

out_sweep <- file.path(tables_dir, "malignant_calling_sweep_runfinal.csv")
write.csv(sweep, out_sweep, row.names = FALSE)

# By-origin summary for two representative cutoffs:
# - strict: all_epi_median_1.5mad (current pipeline)
# - ref_q99: typical paper/forum style using reference 99th percentile
thr_strict <- thresholds[["all_epi_median_1.5mad"]]
thr_ref_q99 <- thresholds[["ref_q99"]]

obs_df <- ann[ann$infercnv_group == "Observation", c("cell_id", "Sample_Origin", "cnv_score")]
obs_df$strict_malignant <- obs_df$cnv_score > thr_strict
obs_df$refq99_malignant <- obs_df$cnv_score > thr_ref_q99

agg <- aggregate(
  obs_df$cnv_score,
  by = list(Sample_Origin = obs_df$Sample_Origin),
  FUN = function(x) c(n = length(x), median = median(x), mean = mean(x))
)
by_origin <- data.frame(
  Sample_Origin = agg$Sample_Origin,
  n = agg$x[, "n"],
  median = agg$x[, "median"],
  mean = agg$x[, "mean"],
  stringsAsFactors = FALSE
)

count_above <- function(flag) {
  t <- aggregate(flag, by = list(Sample_Origin = obs_df$Sample_Origin), FUN = function(x) sum(x, na.rm = TRUE))
  data.frame(Sample_Origin = t$Sample_Origin, count = t$x, stringsAsFactors = FALSE)
}
cs <- count_above(obs_df$strict_malignant); names(cs)[2] <- "strict_malignant"
cr <- count_above(obs_df$refq99_malignant); names(cr)[2] <- "refq99_malignant"
by_origin <- merge(by_origin, cs, by = "Sample_Origin")
by_origin <- merge(by_origin, cr, by = "Sample_Origin")
by_origin$strict_pct <- 100 * by_origin$strict_malignant / by_origin$n
by_origin$refq99_pct <- 100 * by_origin$refq99_malignant / by_origin$n
by_origin <- by_origin[order(-by_origin$n), ]

out_by_origin <- file.path(tables_dir, "malignant_calling_by_origin_runfinal.csv")
write.csv(by_origin, out_by_origin, row.names = FALSE)

# Small base-R boxplot to visualize separation (fast; avoids ggplot deps).
plot_file <- file.path(figures_dir, "cnv_score_boxplot_runfinal.png")
png(plot_file, width = 1400, height = 800)
par(mar = c(9, 5, 4, 2) + 0.1)
lab <- ifelse(ann$infercnv_group == "Ref_Epithelial", "Ref_Epithelial", ann$Sample_Origin)
lab <- factor(lab, levels = c("Ref_Epithelial", sort(unique(ann$Sample_Origin[ann$infercnv_group == "Observation"]))))
boxplot(
  ann$cnv_score ~ lab,
  outline = FALSE,
  las = 2,
  ylab = "CNV score (mean(abs(expr-1)))",
  main = "inferCNV run.final CNV score by group/origin"
)
abline(h = thr_strict, col = "red", lwd = 2)
abline(h = thr_ref_q99, col = "blue", lwd = 2, lty = 2)
legend(
  "topright",
  legend = c(sprintf("strict (all epi median+1.5*MAD)=%.4g", thr_strict), sprintf("ref q99=%.4g", thr_ref_q99)),
  col = c("red", "blue"),
  lwd = 2,
  lty = c(1, 2),
  bty = "n"
)
dev.off()

cat("Saved sweep:", out_sweep, "\n")
cat("Saved by-origin:", out_by_origin, "\n")
cat("Saved plot:", plot_file, "\n")
