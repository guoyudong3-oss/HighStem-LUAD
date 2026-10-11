import sys,importlib.util
sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.metrics import roc_curve, roc_auc_score
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
GENES=["TMSB10","S100A16","HLA-C","S100A14","CKLF","SRPX2"]
spec=importlib.util.spec_from_file_location("f7","07_tcga_genes/38_make_fig7_SIXGENE_O5_copy.py"); f7=importlib.util.module_from_spec(spec); spec.loader.exec_module(f7)
expr=f7.load_expression("data/external/tcga_luad/expression.tsv.gz","data/external/ensembl_to_symbol.csv")
expr_use,merged,gav=f7.prep_data(expr,"data/external/tcga_luad/survival.tsv.gz",GENES)
assert gav==GENES, gav
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
def fp(p): return f"P = {p:.1e}" if p<1e-3 else f"P = {p:.3f}"
def km(fig,x,y,w,h,gene,show_tab=True,ttl=True):
    ax=ax_mm(fig,x,y,w,h); df=merged[[gene,"OS","time_year"]].dropna().copy(); med=df[gene].median(); hi=df[df[gene]>med]; lo=df[df[gene]<=med]
    for g,c,l in [(hi,"#C81D25","High"),(lo,"#128A8A","Low")]:
        k=KaplanMeierFitter().fit(g.time_year,g.OS,label=l); k.plot_survival_function(ax=ax,ci_show=False,color=c,lw=1.0,show_censors=True,censor_styles=dict(marker="|",ms=3,mew=0.5))
    p=logrank_test(hi.time_year,lo.time_year,event_observed_A=hi.OS,event_observed_B=lo.OS).p_value
    ax.set_xlim(0,10); ax.set_ylim(0,1.03); ax.set_xticks([0,2.5,5,7.5,10]); ax.set_xticklabels(["0","2.5","5","7.5","10"]); ax.set_yticks([0,0.5,1.0]); ax.set_yticklabels(["0","0.5","1.0"])
    ax.set_ylabel("Survival probability",fontsize=8,labelpad=1); ax.set_xlabel("Time (years)",fontsize=8,labelpad=1); ax.legend(frameon=False,fontsize=8,loc="upper right",handlelength=1.1,ncol=2,columnspacing=0.6); clean(ax)
    ax.text(0.04,0.06,f"Log-rank p = {p:.3g}",transform=ax.transAxes,fontsize=8)
    if ttl: ax.set_title(gene,fontsize=8,fontweight="bold",pad=2)
    if show_tab:
        tb=ax_mm(fig,x,y+h+11,w,8); tb.set_xlim(0,10); tb.set_ylim(-0.6,1.6); tb.axis("off")
        for yi,(g,l,c) in enumerate([(hi,"High","#C81D25"),(lo,"Low","#128A8A")]):
            for tt in [0,2.5,5,7.5,10]: tb.text(max(tt,0.4),1-yi,str(int((g.time_year>=tt).sum())),fontsize=8,ha="center",va="center")
            tb.text(-0.9,1-yi,l,fontsize=8,ha="right",va="center",color=c)
    return p,len(hi),len(lo)
