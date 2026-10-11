import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
V="results/reproduce_paper/runs/final_malignant_v3_20261009/"
d=pd.read_csv(V+"J1/run_fig5/tables/tcga_highstem_activity_scores_tumor01_fromA.csv")
d=d[np.isfinite(d["OS.time"])&d.OS.notna()].copy(); d["t"]=d["OS.time"]/365; d["OS"]=d.OS.astype(int)
opt=pd.read_csv(V+"J1/run_fig5/tables/Fig5D_summary.csv").cutoff_value[0]; opt=d.PlasticStem_Score.iloc[(d.PlasticStem_Score-opt).abs().argmin()]  # the cutoff in the summary has 15 decimals; snap to the nearest data point to avoid floating-point ties
med=d.PlasticStem_Score.median()
def grp(cut,rule):  # as in the R script: High = score > cutoff
    h=(d.PlasticStem_Score>cut).astype(int); return h
rows=[]
for name,cut in [("optimal cutoff (minprop=0.1)",opt),("median split",med)]:
    h=grp(cut,None); lr=logrank_test(d.t[h==1],d.t[h==0],d.OS[h==1],d.OS[h==0]).p_value
    c=CoxPHFitter().fit(pd.DataFrame({"x":h,"t":d.t,"OS":d.OS}),"t","OS").summary.loc["x"]
    rows.append(dict(method=name,cutoff=cut,n_high=int(h.sum()),n_low=int((1-h).sum()),events_high=int(d.OS[h==1].sum()),events_low=int(d.OS[h==0].sum()),logrank_p=lr,HR_high_vs_low=c["exp(coef)"],HR_CI_low=c["exp(coef) lower 95%"],HR_CI_high=c["exp(coef) upper 95%"],cox_p=c["p"]))
c=CoxPHFitter().fit(d[["PlasticStem_Score","t","OS"]],"t","OS").summary.loc["PlasticStem_Score"]
rows.append(dict(method="continuous Cox (per 1 score unit)",cutoff=np.nan,n_high=np.nan,n_low=np.nan,logrank_p=np.nan,HR_high_vs_low=c["exp(coef)"],HR_CI_low=c["exp(coef) lower 95%"],HR_CI_high=c["exp(coef) upper 95%"],cox_p=c["p"]))
sd=d.PlasticStem_Score.std(); c2=CoxPHFitter().fit(pd.DataFrame({"x":(d.PlasticStem_Score-d.PlasticStem_Score.mean())/sd,"t":d.t,"OS":d.OS}),"t","OS").summary.loc["x"]
rows.append(dict(method="continuous Cox (per 1 SD)",cutoff=np.nan,n_high=np.nan,n_low=np.nan,logrank_p=np.nan,HR_high_vs_low=c2["exp(coef)"],HR_CI_low=c2["exp(coef) lower 95%"],HR_CI_high=c2["exp(coef) upper 95%"],cox_p=c2["p"]))
t=pd.DataFrame(rows); t.insert(1,"n_total",len(d)); t.to_csv(V+"K/K1/Fig5D_three_ways_table.csv",index=False)
pd.set_option("display.width",250); print(t.round(5).to_string(index=False)); print("score SD",round(sd,4))
