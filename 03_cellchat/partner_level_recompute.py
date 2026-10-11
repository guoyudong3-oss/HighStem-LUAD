# Recompute HighStem vs LowStem / DTStem ligand-receptor pair counts from net$count of each CellChat run
import pandas as pd, numpy as np
from scipy.stats import wilcoxon
J="results/reproduce_paper/runs/final_malignant_v3_20261009/J4/"
runs={"main_1200":J+"main_1200/tables/","sensB_fullstem":J+"sensB_fullstem/tables/","sensC_popsizeFALSE":J+"sensC_popsizeFALSE/tables/"}
rows=[];allp=[]
for nm,t in runs.items():
    cm=pd.read_csv(t+"fig3_cellchat_count_matrix.csv",index_col=0)
    partners=[g for g in cm.index if g not in ("HighStem","DTStem","LowStem")]
    def pairs(ref,direction):
        if direction=="outgoing": h,r=cm.loc["HighStem",partners],cm.loc[ref,partners]
        else: h,r=cm.loc[partners,"HighStem"],cm.loc[partners,ref]
        return pd.DataFrame({"run":nm,"partner":partners,"vs":ref,"direction":direction,"high_count":h.values,"ref_count":r.values})
    P={(ref,d):pairs(ref,d) for ref in ("LowStem","DTStem") for d in ("outgoing","incoming")}
    allp+=list(P.values())
    def test(df):
        l2=np.log2((df.high_count+1)/(df.ref_count+1)); 
        try: p=wilcoxon(df.high_count,df.ref_count,alternative="greater").pvalue
        except ValueError: p=np.nan
        return len(df),float(l2.median()),float(2**l2.median()),p,int((df.high_count>df.ref_count).sum())
    combos={"D: outgoing, HighStem vs LowStem":[P[("LowStem","outgoing")]],
            "E: outgoing, HighStem vs DTStem":[P[("DTStem","outgoing")]],
            "pooled D + E(DTStem) outgoing (22 rows)":[P[("LowStem","outgoing")],P[("DTStem","outgoing")]],
            "incoming, HighStem vs LowStem":[P[("LowStem","incoming")]],
            "incoming, HighStem vs DTStem":[P[("DTStem","incoming")]],
            "pooled all four (44 rows, reference only)":list(P.values())}
    for k,v in combos.items():
        n,m,f,p,pos=test(pd.concat(v)); rows.append(dict(run=nm,comparison=k,n_rows=n,median_log2=m,median_fold=f,p_one_sided_signed_rank=p,n_high_gt_ref=pos))
r=pd.DataFrame(rows); r.to_csv(J.replace("J4/","K/K2/")+"partner_level_corrected.csv",index=False)
pd.concat(allp).to_csv(J.replace("J4/","K/K2/")+"partner_pairs_corrected_all.csv",index=False)
pd.set_option("display.width",250); print(r.round(6).to_string(index=False))
