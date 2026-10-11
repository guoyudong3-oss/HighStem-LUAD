
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, rankdata
from scipy.sparse import issparse

# Add parent dir to path for _utils
sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from _utils import setup_logger, ensure_dir

def manual_cytotrace(adata, top_n=200):
    """
    Manually implement CytoTRACE algorithm:
    1. Calculate Gene Counts (GC) per cell.
    2. Correlate every gene expression with GC.
    3. Select top positive correlated genes (Stemness Genes).
    4. Score cells based on mean expression of Stemness Genes.
    5. Smoothing (KNN).
    """
    # 1. Gene Counts
    # Use raw data if available for more accurate counts
    X = adata.raw.X if adata.raw else adata.X
    if issparse(X):
        # binarize
        gene_counts = np.diff(X.tocsr().indptr)
    else:
        gene_counts = np.count_nonzero(X, axis=1)
        
    adata.obs['n_genes_raw'] = gene_counts
    
    # 2. Correlation (Gene vs n_genes)
    # We need dense matrix for correlation usually, or iterate.
    # With 33k cells x 30k genes, dense correlation is heavy.
    # Optimization: Subsample genes to HVG for correlation check? 
    # Or just iterate chunks.
    
    # Let's use HVGs for defining the signature, then project back.
    # This matches CytoTRACE logic (finding variable genes associated with stemness).
    
    # We work on current adata.X (which should be normalized/log1p)
    # If adata.X is sparse:
    if issparse(adata.X):
        data_mat = adata.X.toarray()
    else:
        data_mat = adata.X
        
    # Vectorized correlation
    # corr(X, y) = cov(X, y) / (std(X) * std(y))
    y = gene_counts
    y = y - y.mean()
    y_std = y.std()
    
    correlations = []
    
    # Chunking to save memory if needed, but 128G is enough for 33k*30k float32 (4GB)
    # Let's do it directly.
    X_centered = data_mat - data_mat.mean(axis=0)
    X_std = data_mat.std(axis=0)
    
    # Avoid division by zero
    X_std[X_std == 0] = 1.0 
    
    # Covariance
    # (X.T @ y) / n
    cov = (X_centered.T @ y) / len(y)
    corrs = cov / (X_std * y_std)
    
    gene_cor = pd.Series(corrs, index=adata.var_names)
    
    # 3. Select Stemness Genes
    # CytoTRACE selects genes positively correlated with n_genes
    top_genes = gene_cor.nlargest(top_n).index.tolist()
    
    # 4. Score
    # Mean expression of top genes
    # Geometric mean is used in original CytoTRACE, but arithmetic mean of log-values is equivalent to log of geometric mean.
    # Since data is log1p, mean(log1p) is close.
    # Let's use mean.
    
    score_raw = np.mean(adata[:, top_genes].X, axis=1)
    if issparse(score_raw): score_raw = score_raw.toarray().flatten()
    else: score_raw = np.array(score_raw).flatten()
    
    # 5. Smoothing
    # Use existing neighbors or compute new
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=30)
        
    # Smoothing via connectivity matrix
    # score_smooth = (A @ score) / degree
    # A is weighted adjacency
    conn = adata.obsp['connectivities']
    # Normalize rows to sum to 1
    from sklearn.preprocessing import normalize
    W = normalize(conn, norm='l1', axis=1)
    
    score_smooth = W @ score_raw
    
    # 6. Rank and Normalize (0-1)
    score_final = rankdata(score_smooth) / len(score_smooth)
    
    return score_final, top_genes

def main():
    logger = setup_logger("cytotrace_manual", "results/reproduce_paper/logs/01_c_cytotrace_manual.log")
    
    # Load Data
    input_path = "data/processed/discovery_epithelial_malignant.h5ad"
    logger.info(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    
    # Filter Malignant
    adata_m = adata[adata.obs["malignant_status_recalc"] == "malignant"].copy()
    logger.info(f"Processing {adata_m.n_obs} cells (Manual CytoTRACE)...")
    
    # Run Algorithm
    scores, stem_genes = manual_cytotrace(adata_m)
    logger.info(f"Identified {len(stem_genes)} stemness genes (e.g., {stem_genes[:5]})")
    
    adata_m.obs["cytotrace_score"] = scores
    
    # Define Groups
    q25 = np.quantile(scores, 0.25)
    q75 = np.quantile(scores, 0.75)
    
    def classify(s):
        if s >= q75: return "HighStem"
        if s <= q25: return "LowStem"
        return "DTStem"
        
    adata_m.obs["stem_group"] = [classify(s) for s in scores]
    
    dist = adata_m.obs["stem_group"].value_counts()
    logger.info(f"Groups:\n{dist}")
    
    # Save
    out_path = "data/processed/discovery_malignant_highstem_refined.h5ad"
    adata_m.write(out_path)
    logger.info(f"Saved to {out_path}")

if __name__ == "__main__":
    main()
