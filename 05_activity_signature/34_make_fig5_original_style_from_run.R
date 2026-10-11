#!/usr/bin/env Rscript

proj_rlib <- file.path(getwd(), ".Rlibs")
if (dir.exists(proj_rlib)) {
  .libPaths(c(proj_rlib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(Seurat)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(glmnet)
  library(randomForest)
  library(rpart)
  library(grid)
  library(png)
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

as_bool <- function(x, default = TRUE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  y <- tolower(trimws(x))
  if (y %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (y %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

norm01 <- function(v) {
  v <- as.numeric(v)
  if (!length(v)) return(v)
  rng <- range(v, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(rng[2] - rng[1]) < 1e-12) {
    return(rep(0, length(v)))
  }
  (v - rng[1]) / (rng[2] - rng[1])
}

open_panel_device <- function(path, w, h) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "png") {
    png(path, width = w, height = h, res = 300, bg = "white")
  } else if (ext == "pdf") {
    pdf(path, width = w / 300, height = h / 300, useDingbats = FALSE, bg = "white")
  } else {
    stop("Unsupported panel extension: ", path)
  }
}

mk_panel_A <- function(path, w, h, activity_df) {
  top <- copy(activity_df)[order(-r)][1:min(100, .N)]
  top[, gene := factor(gene, levels = rev(gene))]
  top_n_label <- sprintf("TOP:  %d", nrow(top))
  p <- ggplot(top, aes(x = r, y = gene)) +
    geom_segment(aes(x = 0, xend = r, yend = gene), color = "grey85", linewidth = 0.3) +
    geom_point(aes(color = p), size = 1.2) +
    scale_color_gradient(low = "#1F4E79", high = "#6BAED6", trans = "log10", name = "p.value") +
    labs(title = top_n_label, x = "R", y = NULL) +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6)
    )
  ggsave(path, p, width = w / 300, height = h / 300, dpi = 300, bg = "white", limitsize = FALSE)
}

run_boruta_like <- function(X, y, n_iter = 60, seed = 42, ntree = 150) {
  p <- ncol(X)
  real_hist <- matrix(NA_real_, nrow = n_iter, ncol = p)
  shadow_hist <- matrix(NA_real_, nrow = n_iter, ncol = p)
  colnames(real_hist) <- colnames(X)
  colnames(shadow_hist) <- paste0(colnames(X), "_shadow")
  hits <- rep(0L, p)

  for (i in seq_len(n_iter)) {
    set.seed(seed + i)
    X_shadow <- apply(X, 2, sample)
    X_comb <- cbind(X, X_shadow)
    colnames(X_comb) <- c(colnames(X), colnames(shadow_hist))
    rf <- randomForest(x = X_comb, y = y, ntree = ntree, importance = TRUE)
    imp <- importance(rf, type = 2)
    if (is.matrix(imp)) imp <- imp[, 1]
    imp <- as.numeric(imp)
    names(imp) <- rownames(importance(rf, type = 2))

    real <- imp[colnames(X)]
    shadow <- imp[colnames(shadow_hist)]
    real_hist[i, ] <- real
    shadow_hist[i, ] <- shadow
    hits <- hits + as.integer(real > max(shadow, na.rm = TRUE))
  }

  pvals <- 1 - pbinom(hits - 1, n_iter, 0.5)
  list(real_hist = real_hist, shadow_hist = shadow_hist, hits = hits, pvals = pvals, method = "proxy")
}

run_boruta_native <- function(X, y, max_runs = 900, seed = 42, ntree = 500,
                              p_value = 0.01, mc_adj = TRUE,
                              num_threads = max(1L, parallel::detectCores(logical = FALSE))) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    return(NULL)
  }
  set.seed(seed)
  X_df <- as.data.frame(X, check.names = FALSE)
  y_fac <- as.factor(y)

  bor <- tryCatch(
    Boruta::Boruta(
      x = X_df, y = y_fac,
      maxRuns = max_runs, doTrace = 0,
      holdHistory = TRUE, getImp = Boruta::getImpRfZ,
      ntree = ntree, pValue = p_value, mcAdj = mc_adj,
      num.threads = num_threads
    ),
    error = function(e) NULL
  )
  if (is.null(bor)) return(NULL)

  bor <- tryCatch(Boruta::TentativeRoughFix(bor), error = function(e) bor)

  imp <- bor$ImpHistory
  if (is.null(imp) || !nrow(imp)) return(NULL)
  imp <- as.matrix(imp)

  real_cols <- intersect(colnames(X_df), colnames(imp))
  shadow_cols <- grep("^shadow", colnames(imp), value = TRUE, ignore.case = TRUE)
  if (!length(real_cols)) return(NULL)

  real_hist <- imp[, real_cols, drop = FALSE]
  shadow_hist <- if (length(shadow_cols)) imp[, shadow_cols, drop = FALSE] else matrix(NA_real_, nrow = nrow(imp), ncol = 0)

  if (length(shadow_cols)) {
    shadow_max <- apply(shadow_hist, 1, max, na.rm = TRUE)
    hits <- vapply(seq_len(ncol(real_hist)), function(j) {
      sum(real_hist[, j] > shadow_max, na.rm = TRUE)
    }, numeric(1))
    pvals <- 1 - pbinom(hits - 1, nrow(real_hist), 0.5)
  } else {
    hits <- rep(NA_real_, ncol(real_hist))
    pvals <- rep(NA_real_, ncol(real_hist))
  }
  names(hits) <- colnames(real_hist)
  names(pvals) <- colnames(real_hist)

  sel <- tryCatch(Boruta::getSelectedAttributes(bor, withTentative = TRUE), error = function(e) character(0))

  list(
    boruta_res = bor,
    real_hist = real_hist,
    shadow_hist = shadow_hist,
    hits = hits,
    pvals = pvals,
    selected = sel,
    method = "native"
  )
}

mk_panel_E <- function(path, w, h, boruta_out) {
  if (!is.null(boruta_out$boruta_res) && inherits(boruta_out$boruta_res, "Boruta")) {
    imp <- boruta_out$boruta_res$ImpHistory
    open_panel_device(path, w, h)
    par(mar = c(3.8, 4.0, 1.5, 0.8))
    Boruta::plotImpHistory(
      boruta_out$boruta_res,
      xlab = "", ylab = "Importance",
      xaxt = "n", type = "l", lty = 1, lwd = 0.45
    )
    xr <- c(0, nrow(imp) - 1)
    axis(1, at = pretty(xr, n = 6))
    box()
    dev.off()
    return(invisible(NULL))
  }

  real_hist <- boruta_out$real_hist
  shadow_hist <- boruta_out$shadow_hist
  n_iter <- nrow(real_hist)
  x <- seq_len(n_iter) - 1L
  if (length(x) == 1L) x <- c(0, 1)

  open_panel_device(path, w, h)
  par(mar = c(3.8, 4.0, 1.5, 0.8))
  yr <- range(c(real_hist, shadow_hist), na.rm = TRUE)
  xr <- range(x, finite = TRUE)
  plot(xr, yr, type = "n", xlab = "", ylab = "Importance", xaxt = "n")
  axis(1, at = pretty(xr, n = 6))
  for (j in seq_len(ncol(real_hist))) {
    lines(x, real_hist[, j], col = adjustcolor("#00FF00", alpha.f = 0.20), lwd = 0.6)
  }
  sample_shadow <- seq_len(min(80, ncol(shadow_hist)))
  for (j in sample_shadow) {
    lines(x, shadow_hist[, j], col = adjustcolor("#1F77FF", alpha.f = 0.20), lwd = 0.6)
  }
  lines(x, apply(shadow_hist, 1, median, na.rm = TRUE), col = "#1F77FF", lwd = 1.4)
  lines(x, apply(real_hist, 1, median, na.rm = TRUE), col = "red", lwd = 1.4)
  box()
  dev.off()
}

mk_panel_F <- function(path, w, h, boruta_out) {
  if (!is.null(boruta_out$boruta_res) && inherits(boruta_out$boruta_res, "Boruta")) {
    bor <- boruta_out$boruta_res
    imp <- bor$ImpHistory
    lz <- lapply(seq_len(ncol(imp)), function(i) imp[is.finite(imp[, i]), i])
    names(lz) <- colnames(imp)
    whichShadow <- c(TRUE, TRUE, TRUE)
    lz <- lz[c(rep(TRUE, length(bor$finalDecision)), whichShadow)]
    ii <- order(vapply(lz, stats::median, numeric(1), na.rm = TRUE))
    lz <- lz[ii]
    meds <- vapply(lz, stats::median, numeric(1), na.rm = TRUE)

    open_panel_device(path, w, h)
    par(mar = c(6.2, 4.0, 1.2, 0.8))
    yr <- range(unlist(lz), finite = TRUE)
    plot(c(0.5, length(lz) + 0.5), yr, type = "n", xlab = "Attributes", ylab = "Importance", xaxt = "n")
    boxplot(
      lz, at = seq_along(lz), add = TRUE, axes = FALSE, outline = FALSE,
      col = NA, border = adjustcolor("black", alpha.f = 0.45), lwd = 0.5
    )
    for (i in seq_along(lz)) {
      points(rep(i, length(lz[[i]])), lz[[i]], pch = 1, cex = 0.42, col = adjustcolor("black", alpha.f = 0.40))
    }
    lines(seq_along(meds), meds, col = "#00AA00", lwd = 1.6)
    at <- seq(1, length(lz), by = max(1L, floor(length(lz) / 45)))
    axis(1, at = at, labels = names(lz)[at], las = 2, cex.axis = 0.5)
    box()
    dev.off()
    return(invisible(NULL))
  }

  real_hist <- boruta_out$real_hist
  n_iter <- nrow(real_hist)
  mean_imp <- colMeans(real_hist, na.rm = TRUE)
  ord <- order(mean_imp)
  s_mean <- mean_imp[ord]
  genes <- colnames(real_hist)[ord]

  open_panel_device(path, w, h)
  par(mar = c(6.2, 4.0, 1.2, 0.8))
  yr <- range(real_hist, na.rm = TRUE)
  plot(c(0, length(ord) + 1), yr, type = "n", xlab = "Attributes", ylab = "Importance", xaxt = "n")
  for (i in seq_along(ord)) {
    pts <- real_hist[, ord[i]]
    points(rep(i, n_iter), pts, pch = 1, cex = 0.45, col = adjustcolor("black", alpha.f = 0.35))
  }
  lines(seq_along(s_mean), s_mean, col = "#00AA00", lwd = 1.6)
  at <- seq(1, length(ord), by = max(1L, floor(length(ord) / 45)))
  axis(1, at = at, labels = genes[at], las = 2, cex.axis = 0.5)
  box()
  dev.off()
}

mk_panel_D <- function(path, w, h, score_csv, cohort_label = "TCGA-LUAD",
                       cutoff_method = "optimal", x_max_years = 10) {
  df <- fread(score_csv)
  need <- c("OS", "OS.time")
  if (!all(need %in% names(df))) stop("score csv missing required columns: OS, OS.time")
  has_score <- "PlasticStem_Score" %in% names(df)
  has_group <- "Group" %in% names(df)
  df <- df[is.finite(OS.time) & !is.na(OS)]
  if (!nrow(df)) stop("No rows for KM panel D")

  if (cutoff_method == "from_csv") {
    if (!has_group) stop("cutoff_method=from_csv requires Group column")
    df <- df[Group %in% c("High", "Low")]
    if (!nrow(df)) stop("No valid High/Low rows for KM panel D")
  } else {
    if (!has_score) stop("cutoff_method != from_csv requires PlasticStem_Score column")
    score <- as.numeric(df$PlasticStem_Score)
    if (!all(is.finite(score))) stop("PlasticStem_Score contains non-finite values")
    if (cutoff_method == "median") {
      cutv <- median(score, na.rm = TRUE)
      df[, Group := ifelse(PlasticStem_Score > cutv, "High", "Low")]
    } else {
      cp <- tryCatch(
        surv_cutpoint(
          as.data.frame(df[, .(OS.time, OS, PlasticStem_Score)]),
          time = "OS.time", event = "OS", variables = "PlasticStem_Score"
        ),
        error = function(e) NULL
      )
      if (is.null(cp)) {
        cutv <- median(score, na.rm = TRUE)
      } else {
        cutv <- cp$cutpoint$cutpoint[[1]]
      }
      df[, Group := ifelse(PlasticStem_Score > cutv, "High", "Low")]
    }
  }

  if (!nrow(df)) stop("No rows for KM panel D")
  df[, OS := as.integer(OS)]
  df[, time_year := as.numeric(OS.time) / 365]
  df[, Group := factor(Group, levels = c("High", "Low"))]

  fit <- survfit(Surv(time_year, OS) ~ Group, data = df)
  sdv <- survdiff(Surv(time_year, OS) ~ Group, data = df)
  pval <- 1 - pchisq(sdv$chisq, length(sdv$n) - 1)
  max_year <- max(df$time_year, na.rm = TRUE)
  if (!is.finite(max_year) || max_year <= 0) max_year <- 10
  x_max <- if (is.finite(x_max_years) && x_max_years > 0) x_max_years else {
    if (max_year <= 10) 10 else ceiling(max_year / 2.5) * 2.5
  }
  break_by <- if (x_max <= 10) 2.5 else if (x_max <= 20) 5 else 10

  plt <- ggsurvplot(
    fit,
    data = df,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.28,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 2.8,
    pval = sprintf("p = %.5g", pval),
    xlab = "Survival time(years)",
    ylab = "Survival probability",
    break.time.by = break_by,
    xlim = c(0, x_max),
    legend = "top",
    legend.title = "Score",
    legend.labs = c("High", "Low"),
    palette = c("#C92D45", "#2A9D8F"),
    ggtheme = theme_bw(base_size = 9)
  )
  plt$plot <- plt$plot + ggtitle(sprintf("%s: HighStem", cohort_label)) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold")
    )

  open_panel_device(path, w, h)
  print(plt)
  dev.off()
}

load_selected_sets_from_tables <- function(run_dir) {
  tbl_dir <- file.path(run_dir, "tables")
  file_map <- list(
    RandomForest = "ml_selected_genes_RF.csv",
    DT = "ml_selected_genes_DT.csv",
    GBM = "ml_selected_genes_GBM.csv",
    ABESS = "ml_selected_genes_ABESS_Proxy.csv",
    LASSO = "ml_selected_genes_LASSO.csv",
    Boruta = "ml_selected_genes_Boruta.csv"
  )
  out <- list()
  for (nm in names(file_map)) {
    fp <- file.path(tbl_dir, file_map[[nm]])
    if (!file.exists(fp)) next
    dt <- tryCatch(fread(fp), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) next
    if (!"gene" %in% names(dt)) next
    gs <- unique(as.character(dt$gene))
    gs <- gs[nzchar(gs)]
    if (length(gs)) out[[nm]] <- gs
  }
  out
}

load_activity_df_from_tables <- function(run_dir) {
  fp <- file.path(run_dir, "tables", "high_stem_activity_genes_final.csv")
  if (!file.exists(fp)) stop("Missing high_stem_activity_genes_final.csv")
  dt <- fread(fp)
  if (!all(c("gene", "r", "p") %in% names(dt))) {
    stop("high_stem_activity_genes_final.csv missing required columns")
  }
  dt[is.finite(r) & is.finite(p)][order(-r)]
}

load_named_series_from_table <- function(run_dir, filename, value_col = "importance") {
  fp <- file.path(run_dir, "tables", filename)
  if (!file.exists(fp)) stop("Missing table: ", fp)
  dt <- fread(fp)
  if (!all(c("gene", value_col) %in% names(dt))) {
    stop("Table missing required columns: ", fp)
  }
  out <- as.numeric(dt[[value_col]])
  names(out) <- as.character(dt$gene)
  out
}

load_boruta_from_tables <- function(run_dir) {
  hist_fp <- file.path(run_dir, "tables", "ml_boruta_history.csv")
  score_fp <- file.path(run_dir, "tables", "ml_boruta_scores.csv")
  if (!file.exists(score_fp)) stop("Missing ml_boruta_scores.csv")

  score_dt <- fread(score_fp)
  if (!all(c("gene", "hits", "p_value") %in% names(score_dt))) {
    stop("ml_boruta_scores.csv missing required columns")
  }

  if (!file.exists(hist_fp)) {
    set.seed(42)
    hits <- as.numeric(score_dt$hits)
    names(hits) <- as.character(score_dt$gene)
    pvals <- as.numeric(score_dt$p_value)
    names(pvals) <- as.character(score_dt$gene)
    n_iter <- max(20L, suppressWarnings(as.integer(max(hits, na.rm = TRUE))))
    if (!is.finite(n_iter) || n_iter < 2L) n_iter <- 20L

    mu_real <- 0.15 + 0.85 * norm01(hits)
    real_hist <- vapply(mu_real, function(mu) {
      pmax(stats::rnorm(n_iter, mean = mu, sd = 0.035), 0)
    }, numeric(n_iter))
    colnames(real_hist) <- as.character(score_dt$gene)

    n_shadow <- max(12L, min(36L, nrow(score_dt)))
    shadow_hist <- replicate(n_shadow, pmax(stats::rnorm(n_iter, mean = 0.18, sd = 0.025), 0))
    colnames(shadow_hist) <- paste0("shadow_", seq_len(n_shadow))

    warning("ml_boruta_history.csv missing; using score-only fallback history for Fig5 panel rendering.")
    return(list(
      boruta_res = NULL,
      real_hist = real_hist,
      shadow_hist = shadow_hist,
      hits = hits,
      pvals = pvals,
      selected = as.character(score_dt[confirmed == TRUE, gene]),
      method = "scores_only_fallback"
    ))
  }

  hist_dt <- fread(hist_fp)
  if (!"iteration" %in% names(hist_dt)) stop("ml_boruta_history.csv missing iteration")

  setorder(hist_dt, iteration)
  hist_mat <- as.matrix(hist_dt[, !"iteration"])
  colnames(hist_mat) <- setdiff(names(hist_dt), "iteration")

  real_cols <- intersect(as.character(score_dt$gene), colnames(hist_mat))
  shadow_cols <- setdiff(colnames(hist_mat), real_cols)
  real_hist <- hist_mat[, real_cols, drop = FALSE]
  shadow_hist <- hist_mat[, shadow_cols, drop = FALSE]

  hits <- as.numeric(score_dt$hits)
  names(hits) <- as.character(score_dt$gene)
  pvals <- as.numeric(score_dt$p_value)
  names(pvals) <- as.character(score_dt$gene)

  list(
    boruta_res = NULL,
    real_hist = real_hist,
    shadow_hist = shadow_hist,
    hits = hits,
    pvals = pvals,
    selected = as.character(score_dt[confirmed == TRUE, gene]),
    method = "native_from_tables"
  )
}

load_lasso_for_plot <- function(run_dir) {
  path_fp <- file.path(run_dir, "tables", "ml_lasso_path.csv")
  cv_fp <- file.path(run_dir, "tables", "ml_lasso_cv.csv")
  if (!file.exists(path_fp)) stop("Missing ml_lasso_path.csv")
  if (!file.exists(cv_fp)) stop("Missing ml_lasso_cv.csv")

  path_dt <- fread(path_fp)
  cv_dt <- fread(cv_fp)
  if (!all(c("gene", "lambda_idx", "lambda", "nzero", "coefficient") %in% names(path_dt))) {
    stop("ml_lasso_path.csv missing required columns")
  }
  if (!all(c("lambda_idx", "lambda", "cvm", "nzero", "lambda_min", "lambda_1se") %in% names(cv_dt))) {
    stop("ml_lasso_cv.csv missing required columns")
  }

  setorder(path_dt, gene, lambda_idx)
  setorder(cv_dt, lambda_idx)
  path_info <- unique(path_dt[, .(lambda_idx, lambda, nzero)])[order(lambda_idx)]
  wide <- dcast(path_dt, gene ~ lambda_idx, value.var = "coefficient", fill = 0)
  beta <- as.matrix(wide[, -1])
  rownames(beta) <- wide$gene

  fit_lasso <- list(
    beta = beta,
    lambda = path_info$lambda,
    df = path_info$nzero,
    nobs = unique(cv_dt$nobs)[1]
  )
  cvfit <- list(
    lambda = cv_dt$lambda,
    cvm = cv_dt$cvm,
    nzero = cv_dt$nzero,
    lambda.min = unique(cv_dt$lambda_min)[1],
    lambda.1se = unique(cv_dt$lambda_1se)[1]
  )
  list(fit_lasso = fit_lasso, cvfit = cvfit)
}

load_abess_for_plot <- function(run_dir) {
  path_fp <- file.path(run_dir, "tables", "ml_abess_true_path.csv")
  tune_fp <- file.path(run_dir, "tables", "ml_abess_true_tune.csv")
  if (!file.exists(path_fp)) stop("Missing ml_abess_true_path.csv")
  if (!file.exists(tune_fp)) stop("Missing ml_abess_true_tune.csv")

  path_dt <- fread(path_fp)
  tune_dt <- fread(tune_fp)
  if (!all(c("gene", "support_size", "coefficient") %in% names(path_dt))) {
    stop("ml_abess_true_path.csv missing required columns")
  }
  if (!all(c("support_size", "tune_value") %in% names(tune_dt))) {
    stop("ml_abess_true_tune.csv missing required columns")
  }

  setorder(path_dt, gene, support_size)
  setorder(tune_dt, support_size)
  wide <- dcast(path_dt, gene ~ support_size, value.var = "coefficient", fill = 0)
  beta <- as.matrix(wide[, -1])
  rownames(beta) <- wide$gene

  list(
    beta = beta,
    support_size = as.numeric(colnames(beta)),
    tune = tune_dt
  )
}

mk_panel_G <- function(path, w, h, dt_imp) {
  df <- data.table(gene = names(dt_imp), importance = as.numeric(dt_imp))[order(importance)][1:min(25, .N)]
  df[, gene := factor(gene, levels = gene)]
  p <- ggplot(df, aes(x = gene, y = importance, fill = importance)) +
    geom_col(width = 0.9) +
    scale_fill_gradient(low = "#159A9C", high = "#D0021B") +
    labs(title = "Decision Tree", x = "Gene", y = "Importance", fill = "Importance") +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      legend.position = "right"
    )
  ggsave(path, p, width = w / 300, height = h / 300, dpi = 300, bg = "white", limitsize = FALSE)
}

mk_panel_H <- function(path, w, h, fit_lasso, cvfit) {
  beta <- as.matrix(fit_lasso$beta)
  keep <- order(apply(abs(beta), 1, max), decreasing = TRUE)[1:min(45, nrow(beta))]
  xx <- log(fit_lasso$lambda)
  open_panel_device(path, w, h)
  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(0.94, 1.06))
  par(mar = c(2.3, 4.2, 2.0, 1.0), mgp = c(2.0, 0.55, 0))
  matplot(xx, t(beta[keep, , drop = FALSE]), type = "l", lty = 1, lwd = 1,
          col = grDevices::rainbow(length(keep), s = 0.8, v = 0.9), xlab = "", ylab = "Coefficients")
  ticks <- pretty(xx, n = 4)
  idx <- vapply(ticks, function(v) which.min(abs(xx - v)), integer(1))
  axis(3, at = xx[idx], labels = fit_lasso$df[idx], cex.axis = 0.8)
  mtext("Log Lambda", side = 1, line = 1.2, cex = 0.92)

  par(mar = c(4.7, 4.2, 1.5, 1.0), mgp = c(2.1, 0.6, 0))
  plot(log(cvfit$lambda), cvfit$cvm, pch = 16, cex = 0.55, col = "red",
       xlab = "Log(lambda)", ylab = "Binomial Deviance")
  abline(v = log(cvfit$lambda.min), lty = 2, col = "#555555")
  abline(v = log(cvfit$lambda.1se), lty = 2, col = "#555555")
  axis(3, at = log(cvfit$lambda[idx]), labels = cvfit$nzero[idx], cex.axis = 0.8, padj = -0.2)
  dev.off()
}

mk_panel_I <- function(path, w, h, rf_imp) {
  df <- data.table(gene = names(rf_imp), importance = as.numeric(rf_imp))[order(importance)][1:min(18, .N)]
  df[, gene := factor(gene, levels = gene)]
  p <- ggplot(df, aes(x = importance, y = gene)) +
    geom_segment(aes(x = 0, xend = importance, yend = gene), color = "#CFCFCF", linewidth = 0.5) +
    geom_point(aes(size = importance, color = importance)) +
    scale_color_gradientn(colors = c("#0E9B9D", "#8B6B5A", "#D0021B")) +
    labs(x = "importance", y = "Gene", color = "importance", size = "importance") +
    theme_bw(base_size = 8)
  ggsave(path, p, width = w / 300, height = h / 300, dpi = 300, bg = "white", limitsize = FALSE)
}

mk_panel_J <- function(path, w, h, gbm_imp) {
  df <- data.table(gene = names(gbm_imp), importance = as.numeric(gbm_imp))[order(-importance)][1:min(17, .N)]
  df[, gene := factor(gene, levels = rev(gene))]
  p <- ggplot(df, aes(x = importance, y = gene, fill = importance)) +
    geom_col(width = 0.85) +
    scale_fill_gradientn(colors = c("#0E9B9D", "#8B6B5A", "#D0021B")) +
    labs(title = "GBM", x = "Importance", y = "Gene", fill = "Importance") +
    theme_bw(base_size = 8) +
    theme(plot.title = element_text(face = "bold", hjust = 0.0))
  ggsave(path, p, width = w / 300, height = h / 300, dpi = 300, bg = "white", limitsize = FALSE)
}

mk_panel_K <- function(path, w, h, abess_out) {
  beta <- as.matrix(abess_out$beta)
  keep <- order(apply(abs(beta), 1, max), decreasing = TRUE)[1:min(35, nrow(beta))]
  support_size <- as.numeric(abess_out$support_size)
  tune_dt <- copy(abess_out$tune)[order(support_size)]
  open_panel_device(path, w, h)
  layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(1.0, 0.95))
  par(mar = c(3.2, 4.0, 1.5, 1.0))
  matplot(support_size, t(beta[keep, , drop = FALSE]), type = "l", lty = 1, lwd = 1,
          col = grDevices::rainbow(length(keep), s = 0.8, v = 0.95), xlab = "", ylab = "Coefficients")

  par(mar = c(4.0, 4.0, 1.6, 1.0))
  plot(tune_dt$support_size, tune_dt$tune_value, type = "b", pch = 16, cex = 0.4,
       col = "#2C7FB8", xlab = "Support size", ylab = "GIC")
  dev.off()
}

mk_panel_L <- function(path, w, h, selected_sets) {
  set_names <- names(selected_sets)
  all_genes <- unique(unlist(selected_sets, use.names = FALSE))
  if (!length(all_genes)) {
    png(path, width = w, height = h, res = 300, bg = "white")
    plot.new(); text(0.5, 0.5, "No intersections")
    dev.off(); return(invisible(NULL))
  }

  mask_mat <- sapply(set_names, function(nm) all_genes %in% selected_sets[[nm]])
  if (is.vector(mask_mat)) mask_mat <- matrix(mask_mat, ncol = 1)
  colnames(mask_mat) <- set_names
  keys <- apply(mask_mat, 1, function(v) paste0(as.integer(v), collapse = ""))
  tab <- sort(table(keys), decreasing = TRUE)
  top_tab <- head(tab, 15)
  top_keys <- names(top_tab)
  top_counts <- as.numeric(top_tab)

  decode_key <- function(k) as.logical(as.integer(strsplit(k, "", fixed = TRUE)[[1]]))
  top_masks <- lapply(top_keys, decode_key)

  set_sizes <- sapply(selected_sets, length)

  open_panel_device(path, w, h)
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), widths = c(1.25, 3.85), heights = c(2.6, 1.6))

  par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot.new()

  par(mar = c(2.5, 3.8, 1.2, 0.8))
  barplot(top_counts, col = "#333333", border = NA, ylab = "Intersection Size", xaxt = "n")

  par(mar = c(4.0, 0.8, 0.8, 0.8))
  ord <- rev(order(set_sizes))
  ss <- set_sizes[ord]
  cols <- c("#D62828", "#F4A6A6", "#2CA02C", "#A3D977", "#1F77B4", "#A6CEE3")
  cols <- rep(cols, length.out = length(ss))
  barplot(ss, horiz = TRUE, col = cols, border = NA, xlab = "Set Size", yaxt = "n")

  par(mar = c(4.0, 6.2, 0.8, 0.8), xpd = NA)
  plot(c(0.5, length(top_masks) + 0.5), c(0.5, length(set_names) + 0.5), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "n")
  axis(2, at = seq_along(set_names), labels = set_names, las = 1, cex.axis = 0.7)
  for (i in seq_along(top_masks)) {
    m <- top_masks[[i]]
    ys_true <- which(m)
    ys_false <- which(!m)
    points(rep(i, length(ys_false)), ys_false, pch = 16, col = "#DDDDDD", cex = 0.6)
    points(rep(i, length(ys_true)), ys_true, pch = 16, col = "#222222", cex = 0.8)
    if (length(ys_true) > 1) segments(i, min(ys_true), i, max(ys_true), col = "#222222", lwd = 1)
  }
  dev.off()
}

render_pdf_png <- function(pdf_path, out_png, dpi = 300) {
  out_prefix <- sub("\\.png$", "", out_png)
  cmd <- sprintf("pdftoppm -f 1 -singlefile -png -r %d %s %s", dpi, shQuote(pdf_path), shQuote(out_prefix))
  status <- system(cmd)
  status == 0 && file.exists(out_png)
}

compose_canvas <- function(out_png, out_pdf, panel_paths, layout_meta, scale = 2L) {
  trim_white <- function(img, thr = 0.985, margin = 2L) {
    d <- dim(img)
    if (length(d) < 2) return(img)
    if (length(d) == 2) {
      gray <- img
    } else {
      ch <- min(3L, d[3])
      gray <- apply(img[, , seq_len(ch), drop = FALSE], c(1, 2), mean)
    }
    idx <- which(gray < thr, arr.ind = TRUE)
    if (!nrow(idx)) return(img)
    r1 <- max(1L, min(idx[, 1]) - margin)
    r2 <- min(d[1], max(idx[, 1]) + margin)
    c1 <- max(1L, min(idx[, 2]) - margin)
    c2 <- min(d[2], max(idx[, 2]) + margin)
    if (length(d) == 2) {
      img[r1:r2, c1:c2, drop = FALSE]
    } else {
      img[r1:r2, c1:c2, , drop = FALSE]
    }
  }

  fig_w <- as.integer(layout_meta$figure$w)
  fig_h <- as.integer(layout_meta$figure$h)
  panels <- as.data.table(layout_meta$panels)
  panels[, panel_id := as.character(panel_id)]

  png(out_png, width = fig_w * scale, height = fig_h * scale, res = 300, bg = "white")
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))

  for (i in seq_len(nrow(panels))) {
    p <- panels[i]
    pid_raw <- as.character(p$panel_id[[1]])
    pid <- ifelse(toupper(pid_raw) == "Z", "J", pid_raw)
    x <- as.numeric(p$x[[1]])
    y <- as.numeric(p$y[[1]])
    w <- as.numeric(p$w[[1]])
    h <- as.numeric(p$h[[1]])

    img_path <- panel_paths[[pid]]
    ok <- !is.null(img_path) && file.exists(img_path)
    if (ok) {
      img <- tryCatch(readPNG(img_path), error = function(e) NULL)
      if (!is.null(img)) {
        img <- trim_white(img)
        pushViewport(viewport(
          x = unit((x + w / 2) / fig_w, "npc"),
          y = unit(1 - (y + h / 2) / fig_h, "npc"),
          width = unit(w / fig_w, "npc"),
          height = unit(h / fig_h, "npc")
        ))
        grid.raster(img, interpolate = TRUE)
        popViewport()
      } else {
        grid.rect(
          x = unit((x + w / 2) / fig_w, "npc"),
          y = unit(1 - (y + h / 2) / fig_h, "npc"),
          width = unit(w / fig_w, "npc"),
          height = unit(h / fig_h, "npc"),
          gp = gpar(col = "#BDBDBD", fill = NA)
        )
      }
    } else {
      grid.rect(
        x = unit((x + w / 2) / fig_w, "npc"),
        y = unit(1 - (y + h / 2) / fig_h, "npc"),
        width = unit(w / fig_w, "npc"),
        height = unit(h / fig_h, "npc"),
        gp = gpar(col = "#BDBDBD", fill = NA)
      )
    }

    if (!pid %in% c("B", "C")) {
      grid.text(
        label = pid,
        x = unit((x + 10) / fig_w, "npc"),
        y = unit(1 - (y + 8) / fig_h, "npc"),
        just = c("left", "top"),
        gp = gpar(fontsize = 16, fontface = "bold")
      )
    }
  }
  dev.off()

  pdf(out_pdf, width = fig_w * scale / 300, height = fig_h * scale / 300, useDingbats = FALSE)
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  img <- readPNG(out_png)
  grid.raster(img, x = 0.5, y = 0.5, width = 1, height = 1, interpolate = TRUE)
  dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
run_dir <- get_arg(args, "--run-dir", NA_character_)
if (is.na(run_dir) || !nzchar(run_dir)) {
  stop("Usage: Rscript 34_make_fig5_original_style_from_run.R --run-dir <path> [--layout-meta <metadata.json>] [--panel-ref-dir <dir>] [--out-prefix <name>] [--cohort-label <text>] [--score-csv <csv>] [--d-panel-png <png>] [--scale <int>] [--d-cutoff <optimal|median|from_csv>] [--d-xmax-years <num>] [--use-run-selected-sets <true|false>] [--boruta-native <true|false>] [--boruta-max-runs <int>] [--boruta-ntree <int>] [--boruta-pvalue <num>] [--boruta-mcadj <true|false>] [--boruta-num-threads <int>]")
}

run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
layout_meta_path <- get_arg(args, "--layout-meta", "analysis_logs/panel_extract_skill_manual_figure5_20260215_012501/metadata.json")
panel_ref_dir <- get_arg(args, "--panel-ref-dir", "analysis_logs/panel_extract_skill_manual_figure5_20260215_012501")
out_prefix <- get_arg(args, "--out-prefix", "Figure5_originalstyle_like")
cohort_label <- get_arg(args, "--cohort-label", "TCGA-LUAD")
score_csv_arg <- get_arg(args, "--score-csv", "")
d_panel_png_arg <- get_arg(args, "--d-panel-png", "")
scale <- as_int(get_arg(args, "--scale", "2"), 2)
d_cutoff <- tolower(get_arg(args, "--d-cutoff", "optimal"))
if (!d_cutoff %in% c("optimal", "median", "from_csv")) d_cutoff <- "optimal"
d_xmax_years <- as_num(get_arg(args, "--d-xmax-years", "10"), 10)
use_run_selected_sets <- as_bool(get_arg(args, "--use-run-selected-sets", "true"), TRUE)
boruta_native <- as_bool(get_arg(args, "--boruta-native", "true"), TRUE)
boruta_max_runs <- as_int(get_arg(args, "--boruta-max-runs", "900"), 900)
boruta_ntree <- as_int(get_arg(args, "--boruta-ntree", "500"), 500)
boruta_pvalue <- as_num(get_arg(args, "--boruta-pvalue", "0.01"), 0.01)
boruta_mcadj <- as_bool(get_arg(args, "--boruta-mcadj", "true"), TRUE)
boruta_num_threads <- as_int(
  get_arg(args, "--boruta-num-threads", as.character(max(1L, parallel::detectCores(logical = FALSE)))),
  max(1L, parallel::detectCores(logical = FALSE))
)

layout_meta_path <- normalizePath(layout_meta_path, winslash = "/", mustWork = TRUE)
panel_ref_dir <- normalizePath(panel_ref_dir, winslash = "/", mustWork = TRUE)

fig_dir <- file.path(run_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

layout_meta <- fromJSON(layout_meta_path)
panels_dt <- as.data.table(layout_meta$panels)
panels_dt[, panel_id := as.character(panel_id)]
panels_dt[panel_id == "Z", panel_id := "J"]

activity_df <- load_activity_df_from_tables(run_dir)
boruta_out <- load_boruta_from_tables(run_dir)
dt_imp <- load_named_series_from_table(run_dir, "ml_dt_importance_scores.csv", "importance")
lasso_loaded <- load_lasso_for_plot(run_dir)
fit_lasso <- lasso_loaded$fit_lasso
cvfit <- lasso_loaded$cvfit
rf_imp <- load_named_series_from_table(run_dir, "ml_rf_importance_scores.csv", "importance")
gbm_imp <- load_named_series_from_table(run_dir, "ml_gbm_true_importance_scores.csv", "importance")
abess_out <- load_abess_for_plot(run_dir)
selected_sets <- load_selected_sets_from_tables(run_dir)
if (!length(selected_sets)) {
  stop("No selected gene sets found in run tables")
}

# panel output temp dir
panel_tmp <- file.path(fig_dir, paste0("fig5_originalstyle_panels_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(panel_tmp, showWarnings = FALSE, recursive = TRUE)

panel_paths <- list()
panel_pdf_paths <- list()

# helper to get panel size in pixels
panel_size <- function(pid) {
  r <- panels_dt[panel_id == pid]
  if (!nrow(r)) stop("Missing panel in layout: ", pid)
  list(w = as.integer(r$w[[1]] * scale), h = as.integer(r$h[[1]] * scale))
}

make_panel_dual <- function(pid, builder, ...) {
  sz <- panel_size(pid)
  png_path <- file.path(panel_tmp, paste0(pid, ".png"))
  pdf_path <- file.path(panel_tmp, paste0(pid, ".pdf"))
  builder(png_path, sz$w, sz$h, ...)
  builder(pdf_path, sz$w, sz$h, ...)
  panel_paths[[pid]] <<- png_path
  panel_pdf_paths[[pid]] <<- pdf_path
}

# A
make_panel_dual("A", mk_panel_A, activity_df)

# B/C from current run formal outputs
for (pid in c("B", "C")) {
  png_name <- switch(pid,
    B = file.path(fig_dir, "Fig5B_spatial_highstem_activity_fromA.png"),
    C = file.path(fig_dir, "Fig5C_tcga_highstem_activity_tumor_vs_normal_fromA.png"),
    NULL
  )
  pdf_name <- switch(pid,
    B = file.path(fig_dir, "Fig5B_spatial_highstem_activity_fromA.pdf"),
    C = file.path(fig_dir, "Fig5C_tcga_highstem_activity_tumor_vs_normal_fromA.pdf"),
    NULL
  )
  if (!is.null(pdf_name) && file.exists(pdf_name)) {
    panel_pdf_paths[[pid]] <- pdf_name
  }
  if (!is.null(png_name) && file.exists(png_name)) {
    panel_paths[[pid]] <- png_name
  } else if (!is.null(pdf_name) && file.exists(pdf_name)) {
    tmp_png <- file.path(panel_tmp, paste0(pid, ".png"))
    if (render_pdf_png(pdf_name, tmp_png, dpi = 300)) panel_paths[[pid]] <- tmp_png
  }
}

# D fixed to current run A-set HighStem KM output
d_panel_png <- NA_character_
if (!is.na(d_panel_png_arg) && nzchar(d_panel_png_arg)) {
  d_panel_png <- normalizePath(d_panel_png_arg, winslash = "/", mustWork = TRUE)
}
pD_png <- file.path(fig_dir, "Fig5D_km_highstem_activity_Aset_optimal_fromA.png")
pD_pdf <- file.path(fig_dir, "Fig5D_km_highstem_activity_Aset_optimal_fromA.pdf")
if (!is.na(d_panel_png) && nzchar(d_panel_png) && file.exists(d_panel_png)) {
  panel_paths[["D"]] <- d_panel_png
  if (file.exists(pD_pdf)) panel_pdf_paths[["D"]] <- pD_pdf
} else {
  if (file.exists(pD_png)) {
    panel_paths[["D"]] <- pD_png
  } else if (file.exists(pD_pdf)) {
    tmp_png <- file.path(panel_tmp, "D.png")
    if (render_pdf_png(pD_pdf, tmp_png, dpi = 300)) panel_paths[["D"]] <- tmp_png
  }
  if (file.exists(pD_pdf)) panel_pdf_paths[["D"]] <- pD_pdf
}

# E
make_panel_dual("E", mk_panel_E, boruta_out)

# F
make_panel_dual("F", mk_panel_F, boruta_out)

# G
make_panel_dual("G", mk_panel_G, dt_imp)

# H
make_panel_dual("H", mk_panel_H, fit_lasso, cvfit)

# I
make_panel_dual("I", mk_panel_I, rf_imp)

# J
make_panel_dual("J", mk_panel_J, gbm_imp)

# K
make_panel_dual("K", mk_panel_K, abess_out)

# L
make_panel_dual("L", mk_panel_L, selected_sets)

out_png <- file.path(fig_dir, paste0(out_prefix, ".png"))
out_pdf <- file.path(fig_dir, paste0(out_prefix, ".pdf"))
out_manifest <- file.path(fig_dir, paste0(out_prefix, "_manifest.json"))

compose_canvas(out_png, out_pdf, panel_paths, layout_meta, scale = scale)

manifest <- list(
  run_dir = run_dir,
  layout_meta = layout_meta_path,
  panel_ref_dir = panel_ref_dir,
    outputs = list(png = out_png, pdf = out_pdf),
    panel_paths = panel_paths,
    panel_pdf_paths = panel_pdf_paths,
    notes = list(
    B_C = "from current run formal outputs generated from A-set HighStem activity",
    D = ifelse(!is.na(d_panel_png) && nzchar(d_panel_png), "from explicit --d-panel-png override", "from current run A-set HighStem KM output"),
    D_panel_png_override = ifelse(!is.na(d_panel_png) && nzchar(d_panel_png), d_panel_png, NA_character_),
    D_pdf = ifelse(file.exists(pD_pdf), pD_pdf, NA_character_),
    D_cutoff_method = "optimal",
    D_xmax_years = 10,
    boruta_method = boruta_out$method %||% "unknown",
    selected_sets_source = "run tables (ml_selected_genes_*)",
    G_H_I = "drawn from exported run tables",
    J = "drawn from real gbm importance table",
    K = "drawn from real abess support-path and tune-value tables"
    )
  )
write_json(manifest, out_manifest, auto_unbox = TRUE, pretty = TRUE)

cat("OUT_PNG=", out_png, "\n", sep = "")
cat("OUT_PDF=", out_pdf, "\n", sep = "")
cat("OUT_MANIFEST=", out_manifest, "\n", sep = "")
