# PlasticStem-Lung: Step 2 - Spatial Signature Mapping
# Goal: Calculate HighStem Activity Score (76 genes) for all 56 samples.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(future)
    library(future.apply)
    library(akima)
    library(png)
    library(grid)
    library(jsonlite)
})

# --- Configuration ---
plan("multisession", workers = 10)
options(future.globals.maxSize = 8000 * 1024^2)

OUTPUT_DIR <- "results/reproduce_paper/figures/spatial_module_scores"
DATA_DIR <- "results/reproduce_paper/tables/spatial_scores"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. Load the 76-gene Activity Set
activity_genes <- read.csv("results/reproduce_paper/tables/high_stem_activity_genes_final.csv")$gene
cat("🚀 Loaded", length(activity_genes), "activity genes for scoring.\n")

# 2. Batch Processing Function
calculate_spatial_score <- function(sample_path) {
    sample_name <- basename(sample_path)
    out_pdf <- file.path(OUTPUT_DIR, paste0(sample_name, "_module_score.pdf"))
    out_csv <- file.path(DATA_DIR, paste0(sample_name, "_scores.csv") )
    
    if (file.exists(out_csv)) return(NULL)
    
    tryCatch({
        # A. Load Data
        counts_raw <- Read10X(file.path(sample_path, "filtered_feature_bc_matrix"))
        if (is.list(counts_raw)) {
            counts <- if ("Gene Expression" %in% names(counts_raw)) counts_raw$`Gene Expression` else counts_raw[[1]]
        } else {
            counts <- counts_raw
        }
        
        # B. Load Coords
        pos_path <- file.path(sample_path, "spatial", "tissue_positions.csv")
        if (!file.exists(pos_path)) pos_path <- file.path(sample_path, "spatial", "tissue_positions_list.csv")
        coords <- read.csv(pos_path, header = FALSE)
        if (!grepl("-1", coords[1, 1])) coords <- read.csv(pos_path, header = TRUE)
        colnames(coords)[1:6] <- c("barcode", "tissue", "row", "col", "imagerow", "imagecol")
        rownames(coords) <- coords$barcode
        
        # C. Create Seurat & Score
        common <- intersect(colnames(counts), rownames(coords))
        st_obj <- CreateSeuratObject(counts = counts[, common], assay = "Spatial")
        st_obj <- NormalizeData(st_obj, verbose = FALSE)
        
        # AddModuleScore
        valid_genes <- intersect(activity_genes, rownames(st_obj))
        st_obj <- AddModuleScore(st_obj, features = list(valid_genes), name = "HighStem_Activity_Score")
        
        # D. Save Data
        res_df <- data.frame(
            Barcode = colnames(st_obj),
            Score = st_obj$HighStem_Activity_Score1,
            X_pxl = coords[colnames(st_obj), "imagecol"],
            Y_pxl = coords[colnames(st_obj), "imagerow"],
            In_Tissue = coords[colnames(st_obj), "tissue"]
        )
        write.csv(res_df, out_csv, row.names = FALSE)
        
        # E. Visualization (Smoothed)
        img_path <- file.path(sample_path, "spatial", "tissue_lowres_image.png")
        scale_path <- file.path(sample_path, "spatial", "scalefactors_json.json")
        
        if (file.exists(img_path)) {
            img_data <- readPNG(img_path)
            scales <- fromJSON(scale_path)
            h <- dim(img_data)[1]; w <- dim(img_data)[2]
            sf <- scales$tissue_lowres_scalef
            
            # Prepare Plot DF
            plot_df <- data.frame(
                x = res_df$X_pxl * sf,
                y = h - (res_df$Y_pxl * sf),
                z = res_df$Score
            )
            
            # Interpolation
            grid_x <- seq(0, w, length.out = 300)
            grid_y <- seq(0, h, length.out = 300 * (h/w))
            interp_res <- interp(x = plot_df$x, y = plot_df$y, z = plot_df$z, 
                                 xo = grid_x, yo = grid_y, linear = FALSE, extrap = FALSE)
            df_smooth <- reshape2::melt(interp_res$z, varnames = c("x_idx", "y_idx"), value.name = "val")
            df_smooth$x <- interp_res$x[df_smooth$x_idx]
            df_smooth$y <- interp_res$y[df_smooth$y_idx]
            df_smooth <- na.omit(df_smooth)
            
            g <- rasterGrob(img_data, width = unit(1, "npc"), height = unit(1, "npc"))
            cap <- quantile(df_smooth$val, 0.98, na.rm=T)
            
            p <- ggplot(df_smooth, aes(x = x, y = y)) +
                 annotation_custom(g, xmin = 0, xmax = w, ymin = 0, ymax = h) +
                 geom_raster(aes(fill = val), alpha = 0.6) +
                 scale_fill_gradientn(colors = c("transparent", "yellow", "red"),
                                      values = c(0, 0.3, 1),
                                      limits = c(0, cap), na.value="transparent", oob=scales::squish) +
                 theme_void() +
                 ggtitle(paste0(sample_name, ": Activity Score"))
            
            ggsave(out_pdf, plot = p, width = 8, height = 8)
        }
        
        return(TRUE)
    }, error = function(e) {
        return(FALSE)
    })
}

# 3. Execution
cohort_root <- "data/spatial/GSE307534/Full_Cohort"
sample_dirs <- list.dirs(cohort_root, recursive = TRUE, full.names = TRUE)
valid_samples <- sample_dirs[dir.exists(file.path(sample_dirs, "spatial")) & dir.exists(file.path(sample_dirs, "filtered_feature_bc_matrix"))]

cat("🚀 Calculating Module Scores for 56 samples...\n")
future_lapply(valid_samples, calculate_spatial_score)
cat("✅ All Module Scores Calculated.\n")
