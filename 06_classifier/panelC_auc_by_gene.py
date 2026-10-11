import numpy as np, pandas as pd, anndata as ad
from sklearn.metrics import roc_auc_score
V="results/reproduce_paper/runs/final_malignant_v3_20261009/"
genes=pd.read_csv(V+"M/M1/run/tables/final_consensus_genes.csv").gene.tolist()
a=ad.read_h5ad("data/processed/discovery_malignant.h5ad")
use_raw=a.raw is not None and all(g in set(map(str,a.raw.var_names)) for g in genes)
X=(a.raw[:,genes].X if use_raw else a[:,genes].X); X=X.toarray() if hasattr(X,"toarray") else np.asarray(X)
cnv=pd.read_csv("results/reproduce_paper/runs/paper_main_epiref_latest13k_noclobber_flow_20260213_152245/tables/cnv_paper_exact_cell_metadata_runfinal_ref_mean_3sd.csv").set_index("cell_id")
m13664=set(cnv.index[cnv.malignant_status=="malignant"])
y_new=a.obs_names.astype(str).isin(m13664).astype(int); y_old=(a.obs["malignant_status"].astype(str).values=="malignant").astype(int)
rows=[]
for i,g in enumerate(genes):
    rows.append(dict(gene=g,matrix="adata.raw" if use_raw else "adata.X",
        n_malignant_13664=int(y_new.sum()),n_non_malignant=int((1-y_new).sum()),auc_13664=roc_auc_score(y_new,X[:,i]),
        n_malignant_12902=int(y_old.sum()),n_non_malignant_old=int((1-y_old).sum()),auc_12902_old_label=roc_auc_score(y_old,X[:,i])))
t=pd.DataFrame(rows); t["delta_auc"]=t.auc_13664-t.auc_12902_old_label; t.to_csv(V+"M/M1/panelC_auc_by_gene.csv",index=False)
pd.set_option("display.width",250); print(t.round(4).to_string(index=False))
