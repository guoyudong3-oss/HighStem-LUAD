#!/usr/bin/env Rscript

proj_rlib <- file.path(getwd(), ".Rlibs")
if (dir.exists(proj_rlib)) {
  .libPaths(c(proj_rlib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(glmnet)
  library(randomForest)
  library(rpart)
  library(Boruta)
  library(gbm)
  library(abess)
})

get_arg <- function(args, key, default = NA_character_) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", key, "="), "", hit[[1]]))
  idx <- which(args == key)
  if (length(idx) > 0 && idx[[1]] < length(args)) return(args[idx[[1]] + 1])
  default
}

as_int <- function(x, default) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) default else v
}

as_num <- function(x, default) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) default else v
}

write_gene_list <- function(path, genes) {
  fwrite(data.table(gene = genes), path)
}

copy_if_exists <- function(src, dst_dir) {
  if (file.exists(src)) {
    file.copy(src, file.path(dst_dir, basename(src)), overwrite = TRUE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
source_run <- get_arg(args, "--source-run", "")
out_run <- get_arg(args, "--out-run", "")
if (!nzchar(source_run)) stop("--source-run is required")
if (!nzchar(out_run)) stop("--out-run is required")

source_run <- normalizePath(source_run, winslash = "/", mustWork = TRUE)
if (!dir.exists(out_run)) dir.create(out_run, recursive = TRUE, showWarnings = FALSE)
out_run <- normalizePath(out_run, winslash = "/", mustWork = TRUE)

source_tbl <- file.path(source_run, "tables")
out_tbl <- file.path(out_run, "tables")
out_fig <- file.path(out_run, "figures")
out_log <- file.path(out_run, "logs")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log, recursive = TRUE, showWarnings = FALSE)

boruta_max_runs <- as_int(get_arg(args, "--boruta-max-runs", "900"), 900)
boruta_ntree <- as_int(get_arg(args, "--boruta-ntree", "500"), 500)
boruta_pvalue <- as_num(get_arg(args, "--boruta-pvalue", "0.01"), 0.01)
boruta_num_threads <- as_int(get_arg(args, "--boruta-num-threads", "16"), 16)

needed_files <- c(
  "metadata_strict_latest13k.csv",
  "run_input_provenance.txt",
  "high_stem_activity_genes_final.csv",
  "high_stem_pearson_all_candidates.csv",
  "wgcna_hub_genes.csv",
  "wgcna_hub_genes_candidates_in_adata.csv"
)
for (nm in needed_files) {
  copy_if_exists(file.path(source_tbl, nm), out_tbl)
}

meta_path <- file.path(out_tbl, "metadata_strict_latest13k.csv")
prov_path <- file.path(out_tbl, "run_input_provenance.txt")
activity_path <- file.path(out_tbl, "high_stem_activity_genes_final.csv")
if (!file.exists(meta_path)) stop("Missing metadata_strict_latest13k.csv in source run")
if (!file.exists(prov_path)) stop("Missing run_input_provenance.txt in source run")
if (!file.exists(activity_path)) stop("Missing high_stem_activity_genes_final.csv in source run")

meta <- fread(meta_path)
if (!"cell_id" %in% names(meta)) stop("metadata_strict_latest13k.csv missing cell_id")

prov <- readLines(prov_path, warn = FALSE)
rds_line <- grep("^Input Seurat RDS:", prov, value = TRUE)
if (!length(rds_line)) stop("Cannot parse Input Seurat RDS from run_input_provenance.txt")
rds_path <- sub("^Input Seurat RDS:\\s*", "", rds_line[[1]])
rds_path <- normalizePath(rds_path, winslash = "/", mustWork = TRUE)

message("Loading Seurat object: ", rds_path)
seu <- readRDS(rds_path)
if (!"RNA" %in% names(seu@assays)) stop("RNA assay not found in Seurat object")

expr <- NULL
if ("SCT" %in% names(seu@assays) && "data" %in% Layers(seu[["SCT"]])) {
  DefaultAssay(seu) <- "SCT"
  expr <- GetAssayData(seu, assay = "SCT", layer = "data")
} else if ("data" %in% Layers(seu[["RNA"]])) {
  DefaultAssay(seu) <- "RNA"
  expr <- GetAssayData(seu, assay = "RNA", layer = "data")
} else {
  DefaultAssay(seu) <- "RNA"
  expr <- GetAssayData(seu, assay = "RNA", layer = "counts")
}

activity_df <- fread(activity_path)
if (!all(c("gene", "r", "p") %in% names(activity_df))) {
  stop("high_stem_activity_genes_final.csv missing required columns")
}
activity_genes <- unique(as.character(activity_df$gene))
activity_genes <- activity_genes[activity_genes %in% rownames(expr)]
if (!length(activity_genes)) stop("No activity genes found in expression matrix")

common_cells <- intersect(colnames(seu), meta$cell_id)
meta2 <- copy(meta)[match(common_cells, cell_id)]
train_idx <- meta2$stem_group_strict %in% c("HighStem", "LowStem")
train_cells <- common_cells[train_idx]
train_group <- factor(meta2$stem_group_strict[train_idx], levels = c("LowStem", "HighStem"))
y01 <- as.integer(train_group == "HighStem")

X <- as.matrix(t(expr[activity_genes, train_cells, drop = FALSE]))
X <- scale(X)
X[is.na(X)] <- 0
colnames(X) <- activity_genes

message("Running native Boruta...")
set.seed(42)
bor <- Boruta(
  x = as.data.frame(X, check.names = FALSE),
  y = train_group,
  maxRuns = boruta_max_runs,
  doTrace = 0,
  holdHistory = TRUE,
  getImp = getImpRfZ,
  ntree = boruta_ntree,
  pValue = boruta_pvalue,
  mcAdj = TRUE,
  num.threads = boruta_num_threads
)
bor_fix <- tryCatch(TentativeRoughFix(bor), error = function(e) bor)
bor_selected <- tryCatch(Boruta::getSelectedAttributes(bor_fix, withTentative = TRUE), error = function(e) character(0))
imp <- as.matrix(bor$ImpHistory)
imp_dt <- as.data.table(imp)
imp_dt[, iteration := seq_len(.N)]
setcolorder(imp_dt, c("iteration", setdiff(names(imp_dt), "iteration")))
fwrite(imp_dt, file.path(out_tbl, "ml_boruta_history.csv"))

shadow_cols <- grep("^shadow", colnames(imp), value = TRUE, ignore.case = TRUE)
real_cols <- setdiff(colnames(imp), shadow_cols)
shadow_max <- if (length(shadow_cols)) {
  apply(imp[, shadow_cols, drop = FALSE], 1, max, na.rm = TRUE)
} else {
  rep(NA_real_, nrow(imp))
}
hits <- vapply(real_cols, function(g) sum(imp[, g] > shadow_max, na.rm = TRUE), numeric(1))
pvals <- 1 - pbinom(hits - 1, nrow(imp), 0.5)
final_decision <- rep(NA_character_, length(real_cols))
names(final_decision) <- real_cols
if (!is.null(bor_fix$finalDecision)) {
  fd <- as.character(bor_fix$finalDecision)
  names(fd) <- names(bor_fix$finalDecision)
  final_decision[names(fd)] <- fd
}
boruta_scores_dt <- data.table(
  gene = real_cols,
  hits = as.integer(hits),
  p_value = as.numeric(pvals),
  confirmed = as.logical(pvals < 0.05),
  final_decision = final_decision[real_cols]
)
setorder(boruta_scores_dt, -hits, p_value, gene)
fwrite(boruta_scores_dt, file.path(out_tbl, "ml_boruta_scores.csv"))
write_gene_list(file.path(out_tbl, "ml_selected_genes_Boruta.csv"), sort(unique(bor_selected)))

message("Running native Decision Tree...")
dt_df <- data.frame(y = train_group, X, check.names = FALSE)
set.seed(42)
dt <- rpart(y ~ ., data = dt_df, method = "class", control = rpart.control(cp = 0.001, minsplit = 20))
dt_imp <- dt$variable.importance
if (is.null(dt_imp) || !length(dt_imp)) {
  dt_imp <- setNames(rep(0, ncol(X)), colnames(X))
}
dt_imp_full <- setNames(rep(0, ncol(X)), colnames(X))
dt_imp_full[names(dt_imp)] <- as.numeric(dt_imp)
dt_imp_dt <- data.table(gene = names(dt_imp_full), importance = as.numeric(dt_imp_full))[order(-importance, gene)]
fwrite(dt_imp_dt, file.path(out_tbl, "ml_dt_importance_scores.csv"))
write_gene_list(file.path(out_tbl, "ml_selected_genes_DT.csv"), dt_imp_dt[importance > mean(importance, na.rm = TRUE), gene])

message("Running native LASSO...")
set.seed(42)
fit_lasso <- glmnet(X, y01, family = "binomial", alpha = 1, nlambda = 120, standardize = FALSE)
set.seed(42)
cvfit <- cv.glmnet(X, y01, family = "binomial", alpha = 1, nfolds = 5, standardize = FALSE)
beta <- as.matrix(fit_lasso$beta)
lasso_path_dt <- rbindlist(lapply(seq_along(fit_lasso$lambda), function(i) {
  data.table(
    gene = rownames(beta),
    lambda_idx = i,
    lambda = fit_lasso$lambda[[i]],
    log_lambda = log(fit_lasso$lambda[[i]]),
    nzero = fit_lasso$df[[i]],
    coefficient = beta[, i]
  )
}))
fwrite(lasso_path_dt, file.path(out_tbl, "ml_lasso_path.csv"))
coef_min <- as.matrix(coef(cvfit, s = "lambda.min"))[-1, 1]
lasso_coef_dt <- data.table(
  gene = rownames(as.matrix(coef(cvfit, s = "lambda.min")))[-1],
  coefficient = as.numeric(coef_min)
)
lasso_coef_dt[, abs_coefficient := abs(coefficient)]
setorder(lasso_coef_dt, -abs_coefficient, gene)
fwrite(lasso_coef_dt, file.path(out_tbl, "ml_lasso_coefficients.csv"))
lasso_cv_dt <- data.table(
  lambda_idx = seq_along(cvfit$lambda),
  lambda = cvfit$lambda,
  log_lambda = log(cvfit$lambda),
  cvm = cvfit$cvm,
  nzero = cvfit$nzero,
  lambda_min = cvfit$lambda.min,
  lambda_1se = cvfit$lambda.1se,
  nobs = fit_lasso$nobs
)
fwrite(lasso_cv_dt, file.path(out_tbl, "ml_lasso_cv.csv"))
write_gene_list(file.path(out_tbl, "ml_selected_genes_LASSO.csv"), lasso_coef_dt[abs_coefficient > 1e-5, gene])

message("Running native Random Forest...")
set.seed(42)
rf <- randomForest(x = X, y = train_group, ntree = 300, importance = TRUE)
rf_imp <- importance(rf, type = 2)
if (is.matrix(rf_imp)) rf_imp <- rf_imp[, 1]
rf_imp_dt <- data.table(gene = names(rf_imp), importance = as.numeric(rf_imp))
setorder(rf_imp_dt, -importance, gene)
fwrite(rf_imp_dt, file.path(out_tbl, "ml_rf_importance_scores.csv"))
write_gene_list(file.path(out_tbl, "ml_selected_genes_RF.csv"), rf_imp_dt[importance > mean(importance, na.rm = TRUE), gene])

message("Running native GBM...")
set.seed(42)
gbm_df <- data.frame(y = y01, X, check.names = FALSE)
gbm_fit <- gbm(
  y ~ .,
  data = gbm_df,
  distribution = "bernoulli",
  n.trees = 500,
  interaction.depth = 3,
  shrinkage = 0.01,
  n.minobsinnode = 10,
  verbose = FALSE
)
gbm_imp <- summary(gbm_fit, plotit = FALSE)
gbm_imp_dt <- data.table(gene = as.character(gbm_imp$var), importance = as.numeric(gbm_imp$rel.inf))
setorder(gbm_imp_dt, -importance, gene)
fwrite(gbm_imp_dt, file.path(out_tbl, "ml_gbm_true_importance_scores.csv"))
write_gene_list(file.path(out_tbl, "ml_selected_genes_GBM.csv"), gbm_imp_dt[importance > mean(importance, na.rm = TRUE), gene])
write_gene_list(file.path(out_tbl, "ml_selected_genes_GBM_truth.csv"), gbm_imp_dt[importance > mean(importance, na.rm = TRUE), gene])

message("Running native ABESS...")
set.seed(42)
abess_fit <- abess(X, y01, family = "binomial")
abess_beta <- as.matrix(abess_fit$beta)
abess_support <- as.integer(abess_fit$support.size)
abess_path_dt <- rbindlist(lapply(seq_along(abess_support), function(i) {
  data.table(gene = rownames(abess_beta), support_size = abess_support[[i]], coefficient = abess_beta[, i])
}))
fwrite(abess_path_dt, file.path(out_tbl, "ml_abess_true_path.csv"))
abess_tune_dt <- data.table(
  support_size = abess_support,
  dev = as.numeric(abess_fit$dev),
  tune_value = as.numeric(abess_fit$tune.value),
  edf = as.numeric(abess_fit$edf)
)
abess_tune_dt[, is_best := tune_value == min(tune_value, na.rm = TRUE)]
setorder(abess_tune_dt, support_size)
fwrite(abess_tune_dt, file.path(out_tbl, "ml_abess_true_tune.csv"))
best_idx <- which.min(abess_fit$tune.value)
best_coef <- abess_beta[, best_idx]
abess_coef_dt <- data.table(gene = rownames(abess_beta), coefficient = as.numeric(best_coef))
abess_coef_dt[, abs_coefficient := abs(coefficient)]
setorder(abess_coef_dt, -abs_coefficient, gene)
fwrite(abess_coef_dt, file.path(out_tbl, "ml_abess_true_coefficients.csv"))
abess_selected <- abess_coef_dt[abs_coefficient > 0, gene]
write_gene_list(file.path(out_tbl, "ml_selected_genes_ABESS.csv"), abess_selected)
write_gene_list(file.path(out_tbl, "ml_selected_genes_ABESS_truth.csv"), abess_selected)
write_gene_list(file.path(out_tbl, "ml_selected_genes_ABESS_Proxy.csv"), abess_selected)

selected_sets <- list(
  LASSO = lasso_coef_dt[abs_coefficient > 1e-5, gene],
  Boruta = sort(unique(bor_selected)),
  RF = rf_imp_dt[importance > mean(importance, na.rm = TRUE), gene],
  GBM = gbm_imp_dt[importance > mean(importance, na.rm = TRUE), gene],
  DT = dt_imp_dt[importance > mean(importance, na.rm = TRUE), gene],
  ABESS = abess_selected
)

all_selected <- unlist(selected_sets, use.names = FALSE)
vote_tab <- sort(table(all_selected), decreasing = TRUE)
vote_dt <- data.table(gene = names(vote_tab), votes = as.integer(vote_tab))
fwrite(vote_dt, file.path(out_tbl, "ml_gene_vote_counts.csv"))

cons6 <- vote_dt[votes >= 6, gene]
cons5 <- vote_dt[votes >= 5, gene]
cons4 <- vote_dt[votes >= 4, gene]
final_signature <- if (length(cons6) >= 3) {
  cons6
} else if (length(cons5) >= 3) {
  cons5
} else {
  cons4
}
write_gene_list(file.path(out_tbl, "final_consensus_genes.csv"), final_signature)

summary_dt <- data.table(
  source_run = source_run,
  out_run = out_run,
  activity_gene_count = length(activity_genes),
  train_n = nrow(X),
  high_n = sum(train_group == "HighStem"),
  low_n = sum(train_group == "LowStem"),
  boruta_runs_observed = nrow(imp_dt),
  boruta_max_runs = boruta_max_runs,
  boruta_ntree = boruta_ntree,
  boruta_selected_n = length(selected_sets$Boruta),
  lasso_selected_n = length(selected_sets$LASSO),
  rf_selected_n = length(selected_sets$RF),
  gbm_selected_n = length(selected_sets$GBM),
  dt_selected_n = length(selected_sets$DT),
  abess_selected_n = length(selected_sets$ABESS),
  final_gene_n = length(final_signature)
)
fwrite(summary_dt, file.path(out_tbl, "ml_native_consensus_summary.csv"))

cat("OUT_RUN=", out_run, "\n", sep = "")
cat("FINAL_GENES=", paste(final_signature, collapse = ","), "\n", sep = "")
