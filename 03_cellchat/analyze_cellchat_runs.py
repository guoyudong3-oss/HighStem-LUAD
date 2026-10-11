# Summarise pooled outgoing strength, rank, pair counts and partner-level tests for each CellChat run (one-sided signed-rank, D+E pooled)
import sys, os, numpy as np, pandas as pd
from scipy.stats import wilcoxon
V="results/reproduce_paper/runs/final_malignant_v3_20261009/J4/"
runs={"main_1200":V+"main_1200/tables/","sensB_fullstem":V+"sensB_fullstem/tables/","sensC_popsizeFALSE":V+"sensC_popsizeFALSE/tables/",
      "ref_highratio10":"repro/fig3_prism_light_bundle_20260307/tables/"}
rows=[];pairs=[]
for nm,d in runs.items():
    if not os.path.exists(d+"fig3_cellchat_weight_matrix.csv"): print(nm,"missing or unfinished"); continue
    w=pd.read_csv(d+"fig3_cellchat_weight_matrix.csv",index_col=0)
    out=w.sum(1).sort_values(ascending=False); inn=w.sum(0)
    rk={g:int(list(out.index).index(g))+1 for g in ["HighStem","DTStem","LowStem"]}
    D=pd.read_csv(d+"fig3D_high_vs_low_interaction_pairs.csv"); E=pd.read_csv(d+"fig3E_high_vs_dtstem_interaction_pairs.csv")
    P=pd.concat([D.assign(vs="LowStem"),E.assign(vs="DTStem")])
    l2=np.log2((P.high_count+1)/(P.ref_count+1)); pv=wilcoxon(P.high_count,P.ref_count,alternative="greater").pvalue
    pv_D=wilcoxon(D.high_count,D.ref_count,alternative="greater").pvalue; pv_E=wilcoxon(E.high_count,E.ref_count,alternative="greater").pvalue
    rows.append(dict(run=nm,n_groups=len(w),out_HighStem=out["HighStem"],rank_HighStem=f"{rk['HighStem']}/{len(w)}",out_DTStem=out["DTStem"],rank_DTStem=rk["DTStem"],out_LowStem=out["LowStem"],rank_LowStem=rk["LowStem"],
        in_HighStem=inn["HighStem"],top3_out=";".join(f"{g}:{v:.4f}" for g,v in out.head(3).items()),
        partner_n_rows=len(P),partner_median_log2=float(l2.median()),partner_median_fold=float(2**l2.median()),partner_p_pooled=pv,partner_p_vsLow=pv_D,partner_p_vsDT=pv_E,
        n_high_gt_ref=int((P.high_count>P.ref_count).sum())))
    pairs.append(P.assign(run=nm))
r=pd.DataFrame(rows); r.to_csv(V+"cellchat_runs_summary.csv",index=False); pd.concat(pairs).to_csv(V+"cellchat_partner_pairs_all_runs.csv",index=False)
pd.set_option("display.width",250); print(r.T.to_string())
