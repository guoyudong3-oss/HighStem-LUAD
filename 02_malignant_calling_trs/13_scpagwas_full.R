#!/usr/bin/env Rscript
# ==============================================================================
# scPagwas (full cells) for noclobber run dir
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scPagwas)
  library(data.table)
  library(dplyr)
})

# Patch scPagwas Pathway_block_func to handle empty pathway blocks safely
Pathway_block_func_safe <- function(Pagwas = NULL, Pachrom_block_list, chrom_gwas_list, singlecell = T, celltype = T, chrom_ld, backingpath, n.cores = 1) {
  Pathway_sclm_results <- list()
  Pathway_ld_gwas_data <- list()
  SNP_A <- NULL
  pathway_ckpt_dir <- Sys.getenv("SCPAGWAS_STEP6_DIR", unset = "")
  resume_mode <- tolower(Sys.getenv("SCPAGWAS_RESUME", unset = "1")) %in% c("1", "true", "yes")
  if (nzchar(pathway_ckpt_dir)) {
    dir.create(pathway_ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  }
  max_pathways <- suppressWarnings(as.integer(Sys.getenv("SCPAGWAS_MAX_PATHWAYS", unset = "0")))
  pathway_names <- names(Pachrom_block_list)
  if (!is.na(max_pathways) && max_pathways > 0 && length(pathway_names) > max_pathways) {
    pathway_names <- pathway_names[seq_len(max_pathways)]
    message("Debug mode: limiting pathways to first ", length(pathway_names))
  }
  backingpath_use <- tryCatch(normalizePath(backingpath, winslash = "/", mustWork = FALSE), error = function(e) backingpath)
  if (!dir.exists(backingpath_use)) {
    ok <- tryCatch(dir.create(backingpath_use, recursive = TRUE, showWarnings = FALSE), error = function(e) FALSE)
    if (!isTRUE(ok) || !dir.exists(backingpath_use)) {
      backingpath_use <- file.path(tempdir(), "scPagwas_backing")
      dir.create(backingpath_use, recursive = TRUE, showWarnings = FALSE)
      message("Fallback backingpath: ", backingpath_use)
    }
  }
  # bigstatsr may see NA in this option and fail with "missing value where TRUE/FALSE needed"
  ncores_int <- suppressWarnings(as.integer(n.cores))
  if (is.na(ncores_int) || ncores_int < 1) ncores_int <- 1L
  force_serial <- tolower(Sys.getenv("SCPAGWAS_FORCE_SERIAL", unset = "1")) %in% c("1", "true", "yes")
  if (force_serial && ncores_int != 1L) {
    message("SCPAGWAS_FORCE_SERIAL is ON: forcing n.cores=1 for pathway linking")
    ncores_int <- 1L
    n.cores <- 1L
  }
  current_max <- suppressWarnings(as.integer(getOption("bigstatsr.ncores.max")))
  if (is.na(current_max) || current_max < ncores_int) {
    options(bigstatsr.ncores.max = ncores_int)
  }
  message(paste0("* Start to link gwas and pathway block annotations for ", length(pathway_names), " pathways!"))
  options(bigmemory.allow.dimnames = TRUE)
  pb <- txtProgressBar(style = 3)
  for (pathway in pathway_names) {
    stage <- "init"
    tryCatch({
      ckpt_file <- if (nzchar(pathway_ckpt_dir)) file.path(pathway_ckpt_dir, paste0(pathway, ".rds")) else ""
      if (resume_mode && nzchar(ckpt_file) && file.exists(ckpt_file)) {
        ck <- readRDS(ckpt_file)
        if (singlecell && !is.null(ck$Pathway_sclm_result)) {
          Pathway_sclm_results[[pathway]] <- ck$Pathway_sclm_result
        }
        if (celltype && !is.null(ck$Pathway_ld_gwas_data)) {
          Pathway_ld_gwas_data[[pathway]] <- ck$Pathway_ld_gwas_data
        }
        setTxtProgressBar(pb, which(pathway_names == pathway)/length(pathway_names))
        return(NULL)
      }
      stage <- "load_pathway_blocks"
      Pa_chrom_block <- Pachrom_block_list[[pathway]]
      stage <- "build_chrom_data"
      Pa_chrom_data <- lapply(names(Pa_chrom_block), function(chrom) {
        chrom <- as.character(chrom)
        if (is.na(chrom) || !nzchar(chrom)) {
          return(NULL)
        }
        if (!isTRUE(chrom %in% names(chrom_ld))) {
          return(NULL)
        }
        chrom_block <- Pa_chrom_block[[chrom]]
        ld_data <- chrom_ld[[chrom]]
        if (is.null(chrom_block) || is.null(ld_data) || nrow(chrom_block) == 0 || nrow(ld_data) == 0) {
          return(NULL)
        }
        data.table::setkey(ld_data, SNP_A)
        if (!isTRUE(chrom %in% names(chrom_gwas_list))) {
          warning(paste(chrom, " for gwas is missing, could be a problem!", sep = ""))
          return(NULL)
        }
        if (is.null(chrom_gwas_list[[chrom]])) {
          warning(paste(chrom, " data missing, could be a problem", sep = ""))
          return(NULL)
        }
        rsids <- Pagwas$snp_gene_df[which(Pagwas$snp_gene_df$label %in% chrom_block$label), c("rsid", "label")]
        c2 <- chrom_gwas_list[[chrom]]
        rsids_gwas <- suppressMessages(dplyr::inner_join(rsids, c2))
        if (is.null(nrow(rsids_gwas)) || nrow(rsids_gwas) == 0) {
          return(NULL)
        }
        beta_squared <- rsids_gwas$beta^2
        sub_ld <- ld_data[.(as.vector(rsids_gwas$rsid)), nomatch = 0L]
        return(list(rsids_gwas, beta_squared, sub_ld, chrom_block))
      })
      rm(Pa_chrom_block)
      Pa_chrom_data <- Pa_chrom_data[!sapply(Pa_chrom_data, is.null)]
      if (length(Pa_chrom_data) == 0) {
        setTxtProgressBar(pb, which(pathway_names == pathway)/length(pathway_names))
        return(NULL)
      }
      stage <- "assemble_pa_block"
      pa_block <- list()
      pa_block$block_info <- bigreadr::rbind_df(lapply(seq_len(length(Pa_chrom_data)), function(i) Pa_chrom_data[[i]][[4]]))
      pa_block$snps <- bigreadr::rbind_df(lapply(seq_len(length(Pa_chrom_data)), function(i) Pa_chrom_data[[i]][[1]]))
      pa_block$y <- unlist(lapply(seq_len(length(Pa_chrom_data)), function(i) Pa_chrom_data[[i]][[2]]))
      pa_block$ld_data <- bigreadr::rbind_df(lapply(seq_len(length(Pa_chrom_data)), function(i) Pa_chrom_data[[i]][[3]]))
      rsid_x <- intersect(pa_block$snps$rsid, unique(unlist(pa_block$ld_data[, 1:2])))
      pa_block$ld_data <- as.data.frame(pa_block$ld_data[pa_block$ld_data$SNP_A %in% rsid_x & pa_block$ld_data$SNP_B %in% rsid_x, ])
      rm(Pa_chrom_data)
      if (is.null(pa_block$ld_data) || nrow(pa_block$ld_data) == 0) {
        ld_matrix <- diag(1, nrow = nrow(pa_block$snps))
      }
      else {
        ld_matrix <- scPagwas:::make_ld_matrix(all_snps = pa_block$snps$rsid, ld_data = pa_block$ld_data)
      }
      pa_block$n_snps <- nrow(pa_block$snps)
      pa_block$ld_matrix_squared <- ld_matrix * ld_matrix
      Rn <- paste0(pathway, sample(1:1000, 1))
      cleanup_backing_files <- function(tag) {
        patt <- file.path(backingpath_use, paste0(tag, "*"))
        fs <- Sys.glob(patt)
        if (length(fs) > 0) {
          unlink(fs, recursive = FALSE, force = TRUE)
        }
      }
      cleanup_backing_files(Rn)
      pathway_sclm_result <- NULL
      pathway_ld_result <- NULL
      if (singlecell) {
        stage <- "get_Pathway_sclm"
        pathway_sclm_result <- tryCatch(
          scPagwas:::get_Pathway_sclm(
            pa_block = pa_block,
            pca_scCell_mat = Pagwas$pca_scCell_mat,
            data_mat = Pagwas$data_mat,
            rawPathway_list = Pagwas$rawPathway_list,
            snp_gene_df = Pagwas$snp_gene_df,
            backingpath = backingpath_use,
            n.cores = n.cores,
            Rns = Rn
          ),
          error = function(e) {
            msg <- conditionMessage(e)
            if (ncores_int > 1 && grepl("server socket|cannot be opened|port [0-9]+", msg, ignore.case = TRUE)) {
              message("Pathway ", pathway, ": parallel socket unavailable; retrying with n.cores=1")
              Rn_retry <- paste0(Rn, "_serial")
              cleanup_backing_files(Rn_retry)
              return(scPagwas:::get_Pathway_sclm(
                pa_block = pa_block,
                pca_scCell_mat = Pagwas$pca_scCell_mat,
                data_mat = Pagwas$data_mat,
                rawPathway_list = Pagwas$rawPathway_list,
                snp_gene_df = Pagwas$snp_gene_df,
                backingpath = backingpath_use,
                n.cores = 1,
                Rns = Rn_retry
              ))
            }
            stop(e)
          }
        )
        Pathway_sclm_results[[pathway]] <- pathway_sclm_result
      }
      if (celltype) {
        stage <- "link_pwpca_block"
        pa_block_ct <- tryCatch(
          scPagwas:::link_pwpca_block(pa_block = pa_block, pca_cell_df = data.matrix(Pagwas$pca_cell_df), merge_scexpr = Pagwas$merge_scexpr, snp_gene_df = Pagwas$snp_gene_df, rawPathway_list = Pagwas$rawPathway_list),
          error = function(e) {
            message("Celltype branch skipped for pathway ", pathway, ": ", e$message)
            NULL
          }
        )
        if (!is.null(pa_block_ct)) {
          pathway_ld_result <- pa_block_ct
          Pathway_ld_gwas_data[[pathway]] <- pathway_ld_result
        }
      }
      if (nzchar(ckpt_file)) {
        saveRDS(
          list(
            Pathway_sclm_result = pathway_sclm_result,
            Pathway_ld_gwas_data = pathway_ld_result
          ),
          ckpt_file
        )
        }
      setTxtProgressBar(pb, which(pathway_names == pathway)/length(pathway_names))
    }, error = function(e) {
      call_txt <- tryCatch(paste(deparse(conditionCall(e)), collapse = " "), error = function(ex) "")
      if (!nzchar(call_txt)) {
        message("Skipping pathway ", pathway, " [", stage, "]: ", e$message)
      } else {
        message("Skipping pathway ", pathway, " [", stage, "]: ", e$message, " | call: ", call_txt)
      }
      setTxtProgressBar(pb, which(pathway_names == pathway)/length(pathway_names))
    })
  }
  close(pb)
  if (singlecell) {
    Pathway_sclm_results <- Pathway_sclm_results[!sapply(Pathway_sclm_results, is.null)]
    if (length(Pathway_sclm_results) == 0) {
      stop("No valid pathway single-cell results were generated in step 6.")
    }
    if (length(Pathway_sclm_results) > 0) {
      Pathway_sclm_results <- data.matrix(as.data.frame(Pathway_sclm_results))
      rownames(Pathway_sclm_results) <- colnames(Pagwas$pca_scCell_mat)
      Pagwas$Pathway_sclm_results <- Pathway_sclm_results
    }
  }
  if (celltype) {
    Pagwas$Pathway_ld_gwas_data <- Pathway_ld_gwas_data
  }
  return(Pagwas)
}

use_pathway_patch <- tolower(Sys.getenv("SCPAGWAS_PATCH_PATHWAY", unset = "1")) %in% c("1", "true", "yes")
if (use_pathway_patch) {
  assignInNamespace("Pathway_block_func", Pathway_block_func_safe, ns = "scPagwas")
}

# Patch scPagwas Single_data_input for Seurat v5 compatibility
Single_data_input_safe <- function(Pagwas, Single_data, nfeatures = NULL, Pathway_list = NULL,
                                   assay = "RNA", min_clustercells = 5) {
  message("Input single cell data!")
  if (!("Seurat" %in% class(Single_data))) {
    message("Single_data is not of class Seurat!")
    return(Pagwas)
  }
  Celltype_anno <- data.frame(
    cellnames = rownames(Single_data@meta.data),
    annotation = as.vector(SeuratObject::Idents(Single_data))
  )
  rownames(Celltype_anno) <- Celltype_anno$cellnames
  Afterre_cell_types <- table(Celltype_anno$annotation) > min_clustercells
  Afterre_cell_types <- names(Afterre_cell_types)[Afterre_cell_types]
  message(length(Afterre_cell_types), " cell types are remain, after filter!")
  Celltype_anno <- Celltype_anno[Celltype_anno$annotation %in% Afterre_cell_types, ]
  Single_data <- Single_data[, Celltype_anno$cellnames]
  Pagwas$Celltype_anno <- Celltype_anno

  get_data_layer <- function(obj, assay_name) {
    d1 <- tryCatch(
      SeuratObject::GetAssayData(obj, layer = "data", assay = assay_name),
      error = function(e) e
    )
    if (!inherits(d1, "error")) {
      return(d1)
    }
    d2 <- tryCatch(
      SeuratObject::GetAssayData(obj, slot = "data", assay = assay_name),
      error = function(e) e
    )
    if (!inherits(d2, "error")) {
      return(d2)
    }
    stop("Failed to get assay data layer: ", d1$message, " | ", d2$message)
  }

  data_mat <- NULL
  data_mat <- tryCatch(get_data_layer(Single_data, assay), error = function(e) NULL)
  if (is.null(data_mat) || nrow(data_mat) == 0 || ncol(data_mat) == 0) {
    message("Data layer is empty; running NormalizeData() to create RNA data layer...")
    Single_data <- Seurat::NormalizeData(
      Single_data,
      assay = assay,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
    data_mat <- get_data_layer(Single_data, assay)
  }
  Pagwas$data_mat <- data_mat
  cell_index_by_type <- split(seq_len(nrow(Celltype_anno)), as.character(Celltype_anno$annotation))
  merge_scexpr <- do.call(cbind, lapply(cell_index_by_type, function(idx) {
    Matrix::rowMeans(data_mat[, idx, drop = FALSE])
  }))
  rownames(merge_scexpr) <- rownames(data_mat)
  colnames(merge_scexpr) <- names(cell_index_by_type)
  Pagwas$VariableFeatures <- rownames(Pagwas$data_mat)
  rm(Single_data)
  pagene <- intersect(unique(unlist(Pathway_list)), rownames(Pagwas$data_mat))
  if (length(pagene) < 100) {
    stop("There are little match between rownames of Single_data and pathway genes!")
  }
  Pagwas$VariableFeatures <- intersect(Pagwas$VariableFeatures, pagene)
  if (ncol(merge_scexpr) == 1) {
    a <- colnames(merge_scexpr)
    Pagwas$merge_scexpr <- data.matrix(merge_scexpr[Pagwas$VariableFeatures, ])
    rownames(Pagwas$merge_scexpr) <- Pagwas$VariableFeatures
    colnames(Pagwas$merge_scexpr) <- a
  } else {
    Pagwas$merge_scexpr <- merge_scexpr[Pagwas$VariableFeatures, ]
  }
  return(Pagwas)
}

assignInNamespace("Single_data_input", Single_data_input_safe, ns = "scPagwas")

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
  stop("Usage: Rscript 13_scpagwas_full.R --run-dir <path> [--seurat-rds <path>] [--gwas-file <path>]")
}
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = FALSE)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)

tables_dir <- file.path(run_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
checkpoints_dir <- file.path(run_dir, "checkpoints")
dir.create(checkpoints_dir, recursive = TRUE, showWarnings = FALSE)

# scPagwas expects output.dirs as relative path (it prefixes "./")
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
run_dir_rel <- run_dir
if (startsWith(run_dir, project_root)) {
  run_dir_rel <- sub(paste0("^", project_root, "/?"), "", run_dir)
}

# Paths
input_rds <- get_arg(args, "--seurat-rds", "results/reproduce_paper/r_data/seurat_obj.rds")
gwas_file <- get_arg(args, "--gwas-file", "data/external/gwas/ieu-a-965_full.tsv.gz")
if (!file.exists(input_rds)) stop("Missing Seurat input: ", input_rds)
if (!file.exists(gwas_file)) stop("Missing GWAS file: ", gwas_file)

resume_mode <- tolower(Sys.getenv("SCPAGWAS_RESUME", unset = "1")) %in% c("1", "true", "yes")

step_file <- function(id, name) {
  file.path(checkpoints_dir, sprintf("%02d_%s.rds", as.integer(id), name))
}

restore_or_run <- function(id, name, expr) {
  p <- step_file(id, name)
  if (resume_mode && file.exists(p)) {
    cat("↪ Resume from checkpoint:", p, "\n")
    return(readRDS(p))
  }
  cat("▶ Running step", sprintf("%02d", as.integer(id)), name, "...\n")
  obj <- force(expr)
  saveRDS(obj, p)
  cat("✔ Saved checkpoint:", p, "\n")
  obj
}

cat("================================================================================\n")
cat("🧬 Running scPagwas (Checkpointed)\n")
cat("================================================================================\n")
cat("   input_rds:", normalizePath(input_rds, winslash = "/", mustWork = TRUE), "\n")
cat("   gwas_file:", normalizePath(gwas_file, winslash = "/", mustWork = TRUE), "\n")
cat("   checkpoints:", checkpoints_dir, "\n")
cat("   resume_mode:", ifelse(resume_mode, "ON", "OFF"), "\n")

# 1. Load and Format GWAS Data
gwas_data <- fread(gwas_file)
gwas_formatted <- as.data.frame(gwas_data)
colnames(gwas_formatted) <- tolower(colnames(gwas_formatted))
setnames(
  gwas_formatted,
  old = c("snp", "p", "chr", "bp", "beta", "se"),
  new = c("rsid", "p", "chrom", "pos", "beta", "se"),
  skip_absent = TRUE
)
if (!"maf" %in% colnames(gwas_formatted)) gwas_formatted$maf <- 0.5
gwas_formatted$chrom <- as.character(gwas_formatted$chrom)
gwas_formatted$rsid <- as.character(gwas_formatted$rsid)
gwas_formatted$p <- as.numeric(gwas_formatted$p)
gwas_formatted$pos <- as.numeric(gwas_formatted$pos)
gwas_formatted$beta <- as.numeric(gwas_formatted$beta)
gwas_formatted$se <- as.numeric(gwas_formatted$se)
gwas_formatted$maf <- as.numeric(gwas_formatted$maf)

# Drop rows with missing critical values
keep_cols <- c("rsid", "p", "chrom", "pos", "beta", "se", "maf")
gwas_formatted <- gwas_formatted[complete.cases(gwas_formatted[, keep_cols]), ]
gwas_formatted <- gwas_formatted[gwas_formatted$chrom != "", ]
cat("   GWAS rows after filtering:", nrow(gwas_formatted), "\n")

# 2. Load Seurat (downsample to 15k for stability)
cat("   Loading Seurat object...\n")
seurat_full <- readRDS(input_rds)
seurat_full <- seurat_full[, !is.na(seurat_full$stem_group)]
Idents(seurat_full) <- seurat_full$stem_group
total_cells <- ncol(seurat_full)
cat("   Cells (full):", total_cells, " Genes:", nrow(seurat_full), "\n")

set.seed(42)
target_n <- 15000
target_env <- Sys.getenv("SCPAGWAS_TARGET_N", unset = "")
if (nzchar(target_env)) {
  target_override <- suppressWarnings(as.integer(target_env))
  if (!is.na(target_override)) {
    target_n <- target_override
  }
}
if (target_n > 0 && total_cells > target_n) {
  cells_keep <- sample(colnames(seurat_full), target_n)
  seurat_full <- seurat_full[, cells_keep]
}
cells_used <- ncol(seurat_full)
cat("   Cells (used):", cells_used, "\n")

run_tag <- if (cells_used == total_cells) "full" else paste0("n", cells_used)
output_subdir <- if (run_tag == "full") "scPagwas_full" else "scPagwas"
output_prefix <- if (run_tag == "full") "LUAD_Repro_full" else "LUAD_Repro"
output_dir <- file.path(run_dir_rel, output_subdir)
output_rds <- file.path(tables_dir, paste0("scpagwas_result_", run_tag, ".rds"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

iters_singlecell <- suppressWarnings(as.integer(Sys.getenv("SCPAGWAS_ITERS", unset = "")))
if (is.na(iters_singlecell) || iters_singlecell <= 0) {
  iters_singlecell <- if (run_tag == "full") 5 else 10
}
n_cores <- suppressWarnings(as.integer(Sys.getenv("SCPAGWAS_CORES", unset = "")))
if (is.na(n_cores) || n_cores < 1) {
  n_physical <- suppressWarnings(as.integer(parallel::detectCores(logical = FALSE)))
  n_logical <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
  if (!is.na(n_physical) && n_physical > 0) {
    n_cores <- n_physical
  } else if (!is.na(n_logical) && n_logical > 0) {
    n_cores <- n_logical
  } else {
    n_cores <- 1L
  }
}
force_serial <- tolower(Sys.getenv("SCPAGWAS_FORCE_SERIAL", unset = "1")) %in% c("1", "true", "yes")
if (force_serial) {
  n_cores <- 1L
}
options(bigstatsr.ncores.max = n_cores)
options(mc.cores = n_cores)
for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")) {
  if (force_serial) {
    do.call(Sys.setenv, setNames(list("1"), v))
  } else if (!nzchar(Sys.getenv(v, unset = ""))) {
    do.call(Sys.setenv, setNames(list(as.character(n_cores)), v))
  }
}
cat("   scPagwas output dir:", output_dir, "\n")
cat("   scPagwas prefix:", output_prefix, "\n")
cat("   iters_singlecell:", iters_singlecell, " n.cores:", n_cores, "\n")
cat("   force_serial:", ifelse(force_serial, "ON", "OFF"), "\n")
use_celltype <- tolower(Sys.getenv("SCPAGWAS_CELLTYPE", unset = "0")) %in% c("1", "true", "yes")
cat("   celltype branch:", ifelse(use_celltype, "ON", "OFF"), "\n")

# 3. Load built-in data
data("block_annotation", package = "scPagwas")
data("Genes_by_pathway_kegg", package = "scPagwas")
data("chrom_ld", package = "scPagwas")

Sys.setenv(SCPAGWAS_STEP6_DIR = file.path(checkpoints_dir, "step06_pathways"))
dir.create(Sys.getenv("SCPAGWAS_STEP6_DIR"), recursive = TRUE, showWarnings = FALSE)

cat("   Running scPagwas step-by-step with checkpoints...\n")
tryCatch({
  Pagwas <- restore_or_run(1, "single_data_input", {
    p <- list(rawPathway_list = Genes_by_pathway_kegg)
    p <- Single_data_input_safe(
      Pagwas = p,
      Single_data = seurat_full,
      assay = "RNA",
      Pathway_list = Genes_by_pathway_kegg,
      min_clustercells = 10
    )
    p
  })

  Pagwas <- restore_or_run(2, "pathway_pca", {
    Pathway_pcascore_run(
      Pagwas = Pagwas,
      Pathway_list = Genes_by_pathway_kegg,
      min.pathway.size = 5,
      max.pathway.size = 300
    )
  })

  Pagwas <- restore_or_run(3, "gwas_summary_input", {
    GWAS_summary_input(
      Pagwas = Pagwas,
      gwas_data = gwas_formatted,
      maf_filter = 0.01
    )
  })

  Pagwas <- restore_or_run(4, "snp_to_gene", {
    Pagwas$snp_gene_df <- SnpToGene(
      gwas_data = Pagwas$gwas_data,
      block_annotation = block_annotation,
      marg = 10000
    )
    Pagwas
  })

  Pagwas <- restore_or_run(5, "pathway_annotation", {
    Pathway_annotation_input(Pagwas = Pagwas, block_annotation = block_annotation)
  })

  Pagwas <- restore_or_run(6, "link_pathway_blocks_gwas", {
    Link_pathway_blocks_gwas(
      Pagwas = Pagwas,
      chrom_ld = chrom_ld,
      singlecell = TRUE,
      celltype = use_celltype,
      backingpath = file.path(run_dir, output_subdir, "temp"),
      n.cores = n_cores
    )
  })

  if (use_celltype) {
    Pagwas <- restore_or_run(7, "celltype_contribution", {
      Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
      Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
      Pagwas$bootstrap_results$bp_value <- p.adjust(Pagwas$bootstrap_results$bp_value, method = "fdr")
      utils::write.csv(
        Pagwas$bootstrap_results,
        file = file.path(output_dir, paste0(output_prefix, "_celltypes_bootstrap_results.csv")),
        quote = FALSE
      )
      Pagwas$Pathway_ld_gwas_data <- NULL
      Pagwas$merge_scexpr <- NULL
      Pagwas
    })
  }

  Pagwas <- restore_or_run(8, "singlecell_score", {
    Pagwas$Pathway_ld_gwas_data <- NULL
    scPagwas_perform_score(Pagwas = Pagwas, remove_outlier = TRUE)
  })

  Pagwas <- restore_or_run(9, "gene_pcc", {
    Pagwas$PCC <- scGet_PCC(
      scPagwas.gPAS.score = Pagwas$scPagwas.gPAS.score,
      data_mat = Pagwas$data_mat
    )
    utils::write.csv(
      Pagwas$PCC,
      file = file.path(output_dir, paste0(output_prefix, "_gene_PCC.csv")),
      quote = FALSE
    )
    utils::write.table(
      Pagwas$Pathway_sclm_results,
      file = file.path(output_dir, paste0(output_prefix, "_Pathway_singlecell_lm_results.txt")),
      quote = FALSE,
      sep = "\t"
    )
    # lightweight score export for quick downstream checks
    out_df <- data.frame(
      cell_id = names(Pagwas$scPagwas.gPAS.score),
      scPagwas.gPAS.score = as.numeric(Pagwas$scPagwas.gPAS.score),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      out_df,
      file = file.path(output_dir, paste0(output_prefix, "_singlecell_scPagwas_score_pvalue.Result.csv")),
      row.names = FALSE,
      quote = FALSE
    )
    Pagwas
  })

  saveRDS(Pagwas, output_rds)
  cat("✅ scPagwas checkpointed pipeline complete.\n")
  cat("   Result saved:", output_rds, "\n")
  cat("   PCC saved:", file.path(output_dir, paste0(output_prefix, "_gene_PCC.csv")), "\n")
}, error = function(e) {
  cat("\n❌ Error running scPagwas checkpointed pipeline: ", e$message, "\n")
  traceback(2)
})
