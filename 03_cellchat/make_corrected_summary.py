import pandas as pd, numpy as np
from scipy.stats import wilcoxon
V="results/reproduce_paper/runs/final_malignant_v3_20261009/"
old=pd.read_csv(V+"J4/cellchat_runs_summary.csv")
cor=pd.read_csv(V+"K/K2/partner_level_corrected.csv")
def tst(df):
    l2=np.log2((df.high_count+1)/(df.ref_count+1)); return len(df),float(l2.median()),float(2**l2.median()),float(wilcoxon(df.high_count,df.ref_count,alternative="greater").pvalue),int((df.high_count>df.ref_count).sum())
rows=[]
for _,r in old.iterrows():
    nm=r.run; d=dict(run=nm,n_groups=r.n_groups,out_HighStem=r.out_HighStem,rank_HighStem=r.rank_HighStem,out_DTStem=r.out_DTStem,rank_DTStem=r.rank_DTStem,out_LowStem=r.out_LowStem,rank_LowStem=r.rank_LowStem,top3_out=r.top3_out)
    if nm.startswith("ref_highratio10"):
        B="repro/fig3_prism_light_bundle_20260307/tables/"
        D=pd.read_csv(B+"fig3D_high_vs_low_interaction_pairs.csv"); E=pd.read_csv(B+"fig3E_high_vs_dtstem_interaction_pairs.csv")  # in this bundle the E file compares against DTStem (older script version)
        src="reference bundle (older 18b; E file compares against DTStem)"
    else:
        c=cor[cor.run==nm]; 
        pa=pd.read_csv(V+"K/K2/partner_pairs_corrected_all.csv"); pa=pa[(pa.run==nm)&(pa.direction=="outgoing")]
        D=pa[pa.vs=="LowStem"]; E=pa[pa.vs=="DTStem"]; src="recomputed from net$count (outgoing)"
    for lab,x in (("D_vsLowStem",D),("E_vsDTStem",E)):
        n,m,f,p,pos=tst(x); d.update({f"{lab}_n_partners":n,f"{lab}_median_log2":m,f"{lab}_median_fold":f,f"{lab}_p_one_sided":p,f"{lab}_n_high_gt_ref":pos})
    n,m,f,p,pos=tst(pd.concat([D,E])); d.update(dict(pooled_n_rows=n,pooled_median_log2=m,pooled_median_fold=f,pooled_p_one_sided=p,pooled_note="D and E share the HighStem counts, so rows are not independent; p is descriptive only",source=src))
    rows.append(d)
out=pd.DataFrame(rows); out.to_csv(V+"L/L1/cellchat_runs_summary_corrected.csv",index=False)
pd.set_option("display.width",250); print(out.T.to_string())
