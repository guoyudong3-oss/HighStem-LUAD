#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
})

get_arg <- function(args, key, default = NA_character_) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^", key, "="), "", hit[[1]]))
  idx <- which(args == key)
  if (length(idx) > 0 && idx[[1]] < length(args)) return(args[idx[[1]] + 1])
  default
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

args <- commandArgs(trailingOnly = TRUE)
score_csv <- get_arg(args, "--score-csv", NA_character_)
out_png <- get_arg(args, "--out-png", NA_character_)
out_pdf <- get_arg(args, "--out-pdf", NA_character_)
out_sum <- get_arg(args, "--out-summary", NA_character_)
title_txt <- get_arg(args, "--title", "TCGA-LUAD: Sig3")
cutoff_method <- tolower(get_arg(args, "--cutoff", "from_csv"))
x_max <- as_num(get_arg(args, "--xmax-years", "10"), 10)
show_censor <- as_bool(get_arg(args, "--show-censor", "true"), TRUE)
minprop <- as_num(get_arg(args, "--minprop", "0.1"), 0.1)
if (!is.finite(minprop)) minprop <- 0.1
minprop <- max(0.05, min(0.49, minprop))

if (is.na(score_csv) || !nzchar(score_csv) || !file.exists(score_csv)) {
  stop("Missing --score-csv")
}
if (is.na(out_png) || !nzchar(out_png)) stop("Missing --out-png")
if (is.na(out_pdf) || !nzchar(out_pdf)) stop("Missing --out-pdf")
if (is.na(out_sum) || !nzchar(out_sum)) stop("Missing --out-summary")
if (!cutoff_method %in% c("from_csv", "median", "optimal")) cutoff_method <- "from_csv"

df <- fread(score_csv)
need <- c("OS", "OS.time")
if (!all(need %in% names(df))) stop("score csv missing required columns: OS, OS.time")
if (!"PlasticStem_Score" %in% names(df) && cutoff_method != "from_csv") {
  stop("cutoff method requires PlasticStem_Score")
}

df <- df[is.finite(OS.time) & !is.na(OS)]
df[, OS := as.integer(OS)]
df[, time_year := as.numeric(OS.time) / 365]

if (cutoff_method == "from_csv") {
  if (!"Group" %in% names(df)) stop("cutoff=from_csv requires Group column")
  df <- df[Group %in% c("High", "Low")]
} else if (cutoff_method == "median") {
  cutv <- median(df$PlasticStem_Score, na.rm = TRUE)
  df[, Group := ifelse(PlasticStem_Score > cutv, "High", "Low")]
} else {
  cp <- tryCatch(
    surv_cutpoint(
      as.data.frame(df[, .(OS.time, OS, PlasticStem_Score)]),
      time = "OS.time", event = "OS", variables = "PlasticStem_Score",
      minprop = minprop
    ),
    error = function(e) NULL
  )
  cutv <- if (is.null(cp)) median(df$PlasticStem_Score, na.rm = TRUE) else cp$cutpoint$cutpoint[[1]]
  df[, Group := ifelse(PlasticStem_Score > cutv, "High", "Low")]
}

if (!nrow(df)) stop("No rows after grouping")
df[, Group := factor(Group, levels = c("High", "Low"))]

fit <- survfit(Surv(time_year, OS) ~ Group, data = df)
sdv <- survdiff(Surv(time_year, OS) ~ Group, data = df)
pval <- 1 - pchisq(sdv$chisq, length(sdv$n) - 1)
cox_group <- coxph(Surv(time_year, OS) ~ Group, data = df)
cox_cont <- if ("PlasticStem_Score" %in% names(df)) coxph(Surv(time_year, OS) ~ PlasticStem_Score, data = df) else NULL
cox_group_hr <- unname(exp(coef(cox_group))[["GroupLow"]])
cox_group_p <- unname(summary(cox_group)$coefficients[1, "Pr(>|z|)"])
cox_cont_hr <- if (is.null(cox_cont)) NA_real_ else unname(exp(coef(cox_cont))[["PlasticStem_Score"]])
cox_cont_p <- if (is.null(cox_cont)) NA_real_ else unname(summary(cox_cont)$coefficients[1, "Pr(>|z|)"])

plt <- ggsurvplot(
  fit,
  data = df,
  conf.int = FALSE,
  censor = show_censor,
  censor.shape = 3,
  censor.size = 0.8,
  size = 0.55,
  risk.table = TRUE,
  risk.table.y.text = TRUE,
  risk.table.y.text.col = FALSE,
  risk.table.height = 0.22,
  risk.table.fontsize = 2.3,
  pval = sprintf("p = %.5g", pval),
  pval.size = 3.2,
  pval.coord = c(0.6, 0.16),
  xlab = "Survival time(years)",
  ylab = "Survival probability",
  break.time.by = 2.5,
  xlim = c(0, x_max),
  ylim = c(0, 1.03),
  legend = "top",
  legend.title = "Score",
  legend.labs = c("High", "Low"),
  palette = c("#C8102E", "#1E8E79"),
  ggtheme = theme_classic(base_size = 6.2),
  tables.theme = theme_cleantable(base_size = 6.0)
)

plt$plot <- plt$plot +
  ggtitle(title_txt) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8.2),
    axis.title = element_text(size = 7.0, face = "bold"),
    axis.text = element_text(size = 6.2),
    axis.line = element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(0.08, "cm"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.45),
    legend.title = element_text(face = "bold", size = 7.0),
    legend.text = element_text(size = 6.8),
    legend.key.width = unit(0.45, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    plot.margin = margin(2, 2, 1, 2, unit = "pt")
  )

plt$table <- plt$table +
  theme(
    axis.title.x = element_text(size = 7.0, face = "bold"),
    axis.text.x = element_text(size = 6.2),
    plot.title = element_text(size = 7.8, face = "bold"),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    axis.line = element_blank(),
    plot.margin = margin(0, 2, 1, 2, unit = "pt")
  )

render_scale <- 3
out_w <- 402 * render_scale
out_h <- 421 * render_scale

png(out_png, width = out_w, height = out_h, res = 300, bg = "white", units = "px")
print(plt)
dev.off()

pdf(out_pdf, width = out_w / 300, height = out_h / 300, useDingbats = FALSE)
print(plt)
dev.off()

sum_dt <- data.table(
  signature = if ("signature" %in% names(df)) unique(df$signature)[1] else "custom",
  cutoff_method = cutoff_method,
  minprop = if (cutoff_method == "optimal") minprop else NA_real_,
  cutoff_value = if (exists("cutv")) as.numeric(cutv) else NA_real_,
  N = nrow(df),
  n_high = sum(df$Group == "High", na.rm = TRUE),
  n_low = sum(df$Group == "Low", na.rm = TRUE),
  logrank_p = pval,
  cox_cont_HR = cox_cont_hr,
  cox_cont_p = cox_cont_p,
  cox_group_HR_low_vs_high = cox_group_hr,
  cox_group_p = cox_group_p
)
fwrite(sum_dt, out_sum)

cat("OUT_PNG=", out_png, "\n", sep = "")
cat("OUT_PDF=", out_pdf, "\n", sep = "")
cat("OUT_SUM=", out_sum, "\n", sep = "")
