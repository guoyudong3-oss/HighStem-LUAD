#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
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

must_int <- function(x, name) {
  xi <- suppressWarnings(as.integer(x))
  if (is.na(xi)) stop("Invalid integer for ", name, ": ", x)
  xi
}

must_num <- function(x, name) {
  xn <- suppressWarnings(as.numeric(x))
  if (is.na(xn)) stop("Invalid numeric for ", name, ": ", x)
  xn
}

args <- commandArgs(trailingOnly = TRUE)
run_dir <- get_arg(args, "--run-dir", NULL)
if (is.null(run_dir) || !nzchar(run_dir)) {
  stop("Usage: Rscript 00c_epi_ref_main_step17.R --run-dir <path> [--input-rds <rds>] [--threads <n>] [--cutoff <x>] [--gene-order-file <path>] [--up-to-step <n>]")
}

input_rds <- get_arg(args, "--input-rds", "data/processed/step1_seurat_harmony.rds")
gene_order_arg <- get_arg(args, "--gene-order-file", "")
threads <- must_int(get_arg(args, "--threads", "16"), "--threads")
cutoff <- must_num(get_arg(args, "--cutoff", "0.1"), "--cutoff")
up_to_step <- must_int(get_arg(args, "--up-to-step", "100"), "--up-to-step")
hmm_flag <- tolower(get_arg(args, "--hmm", "true")) %in% c("1", "true", "t", "yes", "y")
denoise_flag <- tolower(get_arg(args, "--denoise", "true")) %in% c("1", "true", "t", "yes", "y")
analysis_mode <- get_arg(args, "--analysis-mode", "subclusters")
no_plot_flag <- tolower(get_arg(args, "--no-plot", "true")) %in% c("1", "true", "t", "yes", "y")
no_prelim_plot_flag <- tolower(get_arg(args, "--no-prelim-plot", "true")) %in% c("1", "true", "t", "yes", "y")
plot_probabilities_flag <- tolower(get_arg(args, "--plot-probabilities", "false")) %in% c("1", "true", "t", "yes", "y")
hmm_report_by <- get_arg(args, "--hmm-report-by", "subcluster")

dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
tables_dir <- file.path(run_dir, "tables")
figures_dir <- file.path(run_dir, "figures")
rds_dir <- file.path(run_dir, "rds")
infercnv_out <- file.path(run_dir, "infercnv_output")
logs_dir <- file.path(run_dir, "logs")
for (d in c(tables_dir, figures_dir, rds_dir, infercnv_out, logs_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(logs_dir, "00c_infercnv_r_main_epi_ref_step17.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(zz)
}, add = TRUE)

cat("============================================================\n")
cat("inferCNV main analysis (epithelial-only reference)\n")
cat("run_dir:", run_dir, "\n")
cat("input_rds:", input_rds, "\n")
cat("threads:", threads, " cutoff:", cutoff, " HMM:", hmm_flag, " denoise:", denoise_flag, " up_to_step:", up_to_step, "\n")
cat("no_plot:", no_plot_flag, " no_prelim_plot:", no_prelim_plot_flag, " plot_probabilities:", plot_probabilities_flag, "\n")
cat("HMM_report_by:", hmm_report_by, "\n")
cat("============================================================\n")
options(scipen = 100)

if (!file.exists(input_rds)) stop("Input RDS not found: ", input_rds)
obj <- readRDS(input_rds)
cat("Loaded Seurat object: cells=", ncol(obj), " genes=", nrow(obj), "\n")

required_cols <- c("Cell_type", "Sample_Origin")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))

cell_type <- as.character(obj$Cell_type)
origin <- as.character(obj$Sample_Origin)

is_epi <- cell_type == "Epithelial cells"
is_normal_origin <- origin %in% c("nLung", "nLN")
ref_epi <- is_epi & is_normal_origin
use_cells <- is_epi

cat("Selected epithelial cells total:", sum(use_cells), "\n")
cat(" - ref epithelial (normal):", sum(ref_epi), "\n")
cat(" - observation epithelial:", sum(is_epi) - sum(ref_epi), "\n")

sub <- obj[, use_cells]
sub_origin <- as.character(sub$Sample_Origin)
sub_is_normal_origin <- sub_origin %in% c("nLung", "nLN")
group <- ifelse(sub_is_normal_origin, "Ref_Epithelial", "Observation")

if (!any(group == "Ref_Epithelial")) stop("No Ref_Epithelial cells found; cannot run epithelial-reference analysis.")

ann <- data.frame(cell_id = colnames(sub), group = group, stringsAsFactors = FALSE)
ann_file <- file.path(tables_dir, "infercnv_annotations_paper_exact.tsv")
write.table(ann[, c("cell_id", "group")], ann_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

gene_order_candidates <- c(
  gene_order_arg,
  "results/reproduce_paper/resources/gencode_v44_gene_pos.txt",
  "data/processed/gencode_v19_gene_pos.txt",
  system.file("extdata", "gencode_v19_gene_pos.txt", package = "infercnv")
)
gene_order_candidates <- gene_order_candidates[nzchar(gene_order_candidates)]
gene_order_file <- NA_character_
for (f in gene_order_candidates) {
  if (file.exists(f)) {
    gene_order_file <- f
    break
  }
}
if (is.na(gene_order_file)) stop("No usable gene order file found.")
cat("Using gene order:", gene_order_file, "\n")

gene_order <- read.table(gene_order_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
if (ncol(gene_order) < 3) stop("Gene order file format invalid: ", gene_order_file)
gene_order <- gene_order[, seq_len(min(4, ncol(gene_order))), drop = FALSE]
gene_order <- gene_order[!duplicated(gene_order[[1]]), , drop = FALSE]
gene_order_dedup <- file.path(tables_dir, "gene_order_dedup_for_infercnv.tsv")
write.table(gene_order, gene_order_dedup, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
gene_order_file <- gene_order_dedup

counts <- tryCatch(
  GetAssayData(sub, layer = "counts", assay = "RNA"),
  error = function(e) GetAssayData(sub, slot = "counts", assay = "RNA")
)
common_genes <- intersect(rownames(counts), gene_order[[1]])
if (length(common_genes) < 200) stop("Too few overlapping genes with gene_order_file: ", length(common_genes))
counts <- counts[common_genes, , drop = FALSE]
cat("Counts after gene-order intersection: genes=", nrow(counts), " cells=", ncol(counts), "\n")

infer_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = ann_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = "Ref_Epithelial"
)

cat("Running infercnv::run ...\n")
infer_obj <- infercnv::run(
  infer_obj,
  cutoff = cutoff,
  out_dir = infercnv_out,
  cluster_by_groups = TRUE,
  denoise = denoise_flag,
  HMM = hmm_flag,
  HMM_report_by = hmm_report_by,
  num_threads = threads,
  analysis_mode = analysis_mode,
  resume_mode = TRUE,
  up_to_step = up_to_step,
  no_plot = no_plot_flag,
  no_prelim_plot = no_prelim_plot_flag,
  plot_probabilities = plot_probabilities_flag
)
cat("infercnv::run done.\n")
cat("Done.\n")
