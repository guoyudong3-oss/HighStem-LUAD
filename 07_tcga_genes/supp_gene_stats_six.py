# Supplementary gene-level statistics: direction, AUC, median-split and continuous Cox models; data loading identical to 38_make_fig7.
import importlib.util, sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
S="07_tcga_genes/38_make_fig7_reference_layout_vector.py"
spec=importlib.util.spec_from_file_location("fig7",S); m=importlib.util.module_from_spec(spec); sys.argv=["x","--run-dir","x","--layout-meta","x"]; spec.loader.exec_module(m)
genes=pd.read_csv("results/reproduce_paper/runs/final_malignant_v3_20261009/J3/fig7/run/tables/final_consensus_genes.csv").gene.tolist()
expr=m.load_expression("data/external/tcga_luad/expression.tsv.gz","data/external/ensembl_to_symbol.csv")
use,merged,ga=m.prep_data(expr,"data/external/tcga_luad/survival.tsv.gz",genes)
rows=[]
for g in ga:
    n=use.loc[use.type=="Normal",g].astype(float); t=use.loc[use.type=="Tumor",g].astype(float)
    y=(use.type=="Tumor").astype(int); auc=roc_auc_score(y,use[g].astype(float))
    d=merged[[g,"OS","time_year"]].dropna().copy()
    hi=d[g]>d[g].median()
    lr=logrank_test(d.loc[hi,"time_year"],d.loc[~hi,"time_year"],event_observed_A=d.loc[hi,"OS"],event_observed_B=d.loc[~hi,"OS"]).p_value
    c=d.rename(columns={g:"x"})[["x","OS","time_year"]]; cph=CoxPHFitter().fit(c,duration_col="time_year",event_col="OS")
    r1=cph.summary.loc["x"]
    c2=pd.DataFrame({"x":hi.astype(int),"OS":d.OS,"time_year":d.time_year}); r2=CoxPHFitter().fit(c2,duration_col="time_year",event_col="OS").summary.loc["x"]
    rows.append(dict(gene=g,n_normal=len(n),n_tumor=len(t),median_normal=n.median(),median_tumor=t.median(),direction="Tumor>Normal" if t.median()>n.median() else "Tumor<Normal",
        mwu_p=mannwhitneyu(n,t).pvalue,auc_tumor_score_raw=auc,n_surv=len(d),logrank_p_median=lr,
        cox_median_HR_high_vs_low=r2["exp(coef)"],cox_median_p=r2["p"],cox_cont_HR_per_log2=r1["exp(coef)"],cox_cont_CI_low=r1["exp(coef) lower 95%"],cox_cont_CI_high=r1["exp(coef) upper 95%"],cox_cont_p=r1["p"]))
df=pd.DataFrame(rows); df.to_csv("results/reproduce_paper/runs/final_malignant_v3_20261009/J3/fig7/gene_level_stats_supplement.csv",index=False)
pd.set_option("display.width",250); print(df.round(5).T.to_string())
