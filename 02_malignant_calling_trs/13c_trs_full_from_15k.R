#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

parse_run_dir <- function(args) {
  run_dir <- NA_character_
  pcc_file <- NA_character_
  input_rds <- "results/reproduce_paper/r_data/seurat_obj.rds"
  out_csv <- NA_character_
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--run-dir" && i < length(args)) {
      run_dir <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--run-dir=")) {
      run_dir <- sub("^--run-dir=", "", a)
      i <- i + 1
      next
    }
    if (a == "--pcc-file" && i < length(args)) {
      pcc_file <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--pcc-file=")) {
      pcc_file <- sub("^--pcc-file=", "", a)
      i <- i + 1
      next
    }
    if (a == "--input-rds" && i < length(args)) {
      input_rds <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--input-rds=")) {
      input_rds <- sub("^--input-rds=", "", a)
      i <- i + 1
      next
    }
    if (a == "--out-csv" && i < length(args)) {
      out_csv <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (startsWith(a, "--out-csv=")) {
      out_csv <- sub("^--out-csv=", "", a)
      i <- i + 1
      next
    }
    i <- i + 1
  }
  if (is.na(run_dir) || !nzchar(run_dir)) {
    stop(
      "Usage: Rscript 13c_trs_full_from_15k.R --run-dir <path> ",
      "[--pcc-file <csv>] [--input-rds <rds>] [--out-csv <csv>]"
    )
  }
  list(
    run_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    pcc_file = if (!is.na(pcc_file) && nzchar(pcc_file)) normalizePath(pcc_file, winslash = "/", mustWork = FALSE) else NA_character_,
    input_rds = normalizePath(input_rds, winslash = "/", mustWork = FALSE),
    out_csv = if (!is.na(out_csv) && nzchar(out_csv)) normalizePath(out_csv, winslash = "/", mustWork = FALSE) else NA_character_
  )
}

cli <- parse_run_dir(commandArgs(trailingOnly = TRUE))
run_dir <- cli$run_dir
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

tables_dir <- file.path(run_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

pcc_candidates <- unique(na.omit(c(
  cli$pcc_file,
  file.path(run_dir, "scPagwas", "LUAD_Repro_gene_PCC.csv"),
  file.path(run_dir, "scPagwas_full", "LUAD_Repro_full_gene_PCC.csv")
)))
pcc_path <- pcc_candidates[file.exists(pcc_candidates)][1]
if (is.na(pcc_path) || !nzchar(pcc_path)) {
  stop("Missing PCC file. Tried: ", paste(pcc_candidates, collapse = ", "))
}

input_rds <- cli$input_rds
if (!file.exists(input_rds)) {
  stop(paste0("Missing Seurat object: ", input_rds))
}

top_n <- suppressWarnings(as.integer(Sys.getenv("SCPAGWAS_TOP_N", unset = "1000")))
if (is.na(top_n) || top_n <= 0) {
  top_n <- 1000
}

output_csv <- if (!is.na(cli$out_csv) && nzchar(cli$out_csv)) cli$out_csv else file.path(tables_dir, "scpagwas_trs_scores_full_from_15k.csv")

cat("================================================================================\n")
cat("🧬 TRS module score (full coverage from 15k scPagwas genes)\n")
cat("================================================================================\n")

cat("   Loading PCC file...\n")
pcc <- read.csv(pcc_path, stringsAsFactors = FALSE)
if (!all(c("X", "PCC") %in% colnames(pcc))) {
  stop("PCC file missing required columns X and PCC")
}

pcc_desc <- pcc[order(pcc$PCC, decreasing = TRUE), ]
top_genes <- pcc_desc$X[seq_len(min(top_n, nrow(pcc_desc)))]
pcc_asc <- pcc[order(pcc$PCC, decreasing = FALSE), ]
down_genes <- pcc_asc$X[seq_len(min(top_n, nrow(pcc_asc)))]

cat("   Loading Seurat object...\n")
seu <- readRDS(input_rds)
seu <- seu[, !is.na(seu$stem_group)]
cat("   Cells:", ncol(seu), " Genes:", nrow(seu), "\n")

gene_universe <- rownames(seu)
top_genes <- intersect(top_genes, gene_universe)
down_genes <- intersect(down_genes, gene_universe)
cat("   Top genes in data:", length(top_genes), "\n")
cat("   Down genes in data:", length(down_genes), "\n")

cat("   Computing module scores...\n")
seu <- AddModuleScore(seu, features = list(top_genes), name = "scPagwas.TRS.Score")
seu <- AddModuleScore(seu, features = list(down_genes), name = "scPagwas.downTRS.Score")

df <- data.frame(
  cell_id = colnames(seu),
  trs_score = seu$scPagwas.TRS.Score1,
  down_trs_score = seu$scPagwas.downTRS.Score1,
  stringsAsFactors = FALSE
)
write.csv(df, output_csv, row.names = FALSE)
cat("✅ TRS module scores saved to:", output_csv, "\n")

trs_15k <- file.path(tables_dir, "scpagwas_trs_scores.csv")
if (file.exists(trs_15k)) {
  cat("   Comparing to 15k TRS scores...\n")
  df15 <- read.csv(trs_15k, stringsAsFactors = FALSE)
  merged <- merge(df, df15, by = "cell_id")
  if (nrow(merged) > 5) {
    r <- suppressWarnings(cor(merged$trs_score.x, merged$trs_score.y))
    cat("   Correlation (full module score vs 15k TRS):", r, "\n")
  }
}
