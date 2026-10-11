# Fig S5: Boruta, decision tree, random forest, GBM, LASSO and ABESS proxy
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/FigS5/"; T=J+"P/P2/run_fig5/tables/"
fig=new_fig(222)
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
# Boruta importance history
h=pd.read_csv(T+"ml_boruta_history.csv"); a=ax_mm(fig,14,8,70,50); gl=[c for c in h.columns if c!="iteration"]
for g in gl: a.plot(h.iteration,h[g],color="#9BE08A",lw=0.4)
a.set_xlabel("Boruta iteration",fontsize=8,labelpad=1); a.set_title("Boruta: importance per round",fontsize=8,fontweight="bold",pad=2); a.set_ylabel("Importance",fontsize=8,labelpad=1); a.set_xlim(1,h.iteration.max()); clean(a)
# Boruta hits per gene
bs=pd.read_csv(T+"ml_boruta_scores.csv").sort_values("hits",ascending=False).reset_index(drop=True); b=ax_mm(fig,106,8,70,50)
b.bar(range(len(bs)),bs.hits,color=np.where(bs.confirmed,"#3C9A5F","#BBBBBB"),width=0.9); b.set_xlabel("Genes (ranked)",fontsize=8,labelpad=1); b.set_ylabel("Boruta hits (of 20 rounds)",fontsize=8,labelpad=1); b.set_xticks([]); clean(b); b.set_title("Boruta: hits per gene (green, confirmed)",fontsize=8,fontweight="bold",pad=2)
# decision-tree importance (top genes)
def topbar(ax,f,col,n,color,title):
    d=pd.read_csv(T+f).sort_values("score",ascending=False).head(n).iloc[::-1]
    ax.barh(range(n),d.score,color=color,height=0.8); ax.set_yticks(range(n)); ax.set_yticklabels(d.gene,fontsize=8); ax.tick_params(axis="y",length=0,pad=1.5)
    ax.set_xlabel("Importance",fontsize=8,labelpad=1); ax.set_title(title,fontsize=8,fontweight="bold",pad=2); clean(ax)
c=ax_mm(fig,16,76,36,66); topbar(c,"ml_DT_scores.csv","score",20,"#2A9D9D","Decision tree (top 20)")
d_=ax_mm(fig,76,76,36,66); topbar(d_,"ml_RF_scores.csv","score",20,"#C0504D","Random forest (top 20)")
e=ax_mm(fig,136,76,36,66); topbar(e,"ml_GBM_scores.csv","score",20,"#8C6D31","GBM (top 20)")
# LASSO path and CV
lp=pd.read_csv(T+"ml_lasso_path.csv"); cv=pd.read_csv(T+"ml_lasso_cv.csv"); f1=ax_mm(fig,14,162,70,20); f2=ax_mm(fig,14,190,70,20)
for g,m in lp.groupby("gene"): f1.plot(np.log(m["lambda"]),m.coefficient,lw=0.4)
f1.set_ylabel("Coefficient",fontsize=8,labelpad=1); f1.set_xticks([]); clean(f1)
f2.errorbar(np.log(cv["lambda"]),cv.cvm,yerr=cv.cvsd,fmt="o",ms=1.5,color="red",ecolor="#F2B8B8",lw=0.4)
for lam in (cv.lambda_min.iloc[0],cv.lambda_1se.iloc[0]): f2.axvline(np.log(lam),color="#555555",ls="--",lw=0.6)
f2.set_xlabel("log(lambda)",fontsize=8,labelpad=1); f1.set_title("LASSO",fontsize=8,fontweight="bold",pad=2); f2.set_ylabel("CV MSE",fontsize=8,labelpad=1); clean(f2)
# ABESS proxy (RFE-SVC)
ap=pd.read_csv(T+"ml_abess_proxy_path.csv"); tu=pd.read_csv(T+"ml_abess_proxy_tune.csv").sort_values("support_size"); g1=ax_mm(fig,106,162,70,20); g2=ax_mm(fig,106,190,70,20)
for g,m in ap.groupby("gene"): g1.plot(m.support_size,m.coefficient,lw=0.4)
g1.set_ylabel("Coefficient",fontsize=8,labelpad=1); g1.set_xticks([]); clean(g1)
g2.plot(tu.support_size,tu.tune_value,"-o",ms=1.5,lw=0.7,color="#2C7FB8"); g2.set_xlabel("Number of genes retained",fontsize=8,labelpad=1); g1.set_title("ABESS proxy (RFE with linear SVM)",fontsize=8,fontweight="bold",pad=2); g2.set_ylabel("Training error",fontsize=8,labelpad=1); clean(g2)
for t,x,y in [("a",1,1),("b",1,70),("c",61,70),("d",121,70),("e",1,155),("f",93,155)]: letter(fig,x,y,t)
out=V+"FigS5_v2"; save(fig,out); print(check(out+".pdf",V+"FigS5_check.json"))
