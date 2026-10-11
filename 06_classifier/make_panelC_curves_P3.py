# Export single-gene ROC curve points (13,664 inferCNV malignant labels, adata.raw, six genes) for the single-gene ROC panel
import numpy as np, pandas as pd, anndata as ad
from sklearn.metrics import roc_curve, auc
V="results/reproduce_paper/runs/final_malignant_v3_20261009/"
genes=pd.read_csv(V+"M/M1/run/tables/final_consensus_genes.csv").gene.tolist()
a=ad.read_h5ad("data/processed/discovery_malignant.h5ad")
X=a.raw[:,genes].X; X=X.toarray() if hasattr(X,"toarray") else np.asarray(X)
cnv=pd.read_csv("results/reproduce_paper/runs/paper_main_epiref_latest13k_noclobber_flow_20260213_152245/tables/cnv_paper_exact_cell_metadata_runfinal_ref_mean_3sd.csv").set_index("cell_id")
y=a.obs_names.astype(str).isin(set(cnv.index[cnv.malignant_status=="malignant"])).astype(int)
out={}
for i,g in enumerate(genes):
    f,t,_=roc_curve(y,X[:,i]); out[g+"_fpr"]=f; out[g+"_tpr"]=t; print(g,round(auc(f,t),3))
np.savez(V+"P/P3/panelC_roc_curves.npz",genes=np.array(genes),**out)
