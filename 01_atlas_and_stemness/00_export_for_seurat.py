
import os
import scanpy as sc
import pandas as pd
import scipy.io
from scipy import sparse
import shutil

def main():
    print("Loading h5ad...")
    # Update to the REFINED malignant object
    adata = sc.read_h5ad("data/processed/discovery_malignant_highstem_refined.h5ad")
    
    out_dir = "results/reproduce_paper/r_data"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    print(f"Exporting to {out_dir} for Seurat...")
    
    # 1. Counts Matrix (genes x cells)
    # Ensure we use raw counts for WGCNA
    if adata.raw:
        data = adata.raw.X
        var_names = adata.raw.var_names
    else:
        data = adata.X
        var_names = adata.var_names
        
    if sparse.issparse(data):
        data = data.T # Scanpy is cells x genes, Seurat wants genes x cells usually? 
        # Actually Matrix::readMM reads it as is. 
        # But 10x convention: matrix.mtx is (features x barcodes).
        # Scanpy X is (barcodes x features). So we must transpose.
        scipy.io.mmwrite(os.path.join(out_dir, "matrix.mtx"), data)
    else:
        scipy.io.mmwrite(os.path.join(out_dir, "matrix.mtx"), sparse.csr_matrix(data).T)
        
    # 2. Barcodes (cells)
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as f:
        for bc in adata.obs_names:
            f.write(f"{bc}\n")
            
    # 3. Features (genes)
    # Seurat Read10X expects: id \t name \t type (optional)
    # We just write names
    pd.DataFrame(var_names).to_csv(os.path.join(out_dir, "features.tsv"), sep="\t", header=False, index=False)
    
    # 4. Metadata (obs)
    adata.obs.to_csv(os.path.join(out_dir, "metadata.csv"))
    
    print("Export done.")

if __name__ == "__main__":
    main()
