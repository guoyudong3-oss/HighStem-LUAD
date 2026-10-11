# Fig S8: six genes vs PDPK1/AKT1/MTOR scatter plots (Spearman; n = 517)
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
from scipy.stats import spearmanr
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/FigS9/"
d=pd.read_csv(J+"L/L2/cache/fig8_bioinfo_tcga_tumor_log2tpm_core_pi3k.csv"); d.columns=[c.replace("HLA.C","HLA-C") for c in d.columns]; assert len(d)==517
genes=["TMSB10","S100A16","HLA-C","S100A14","CKLF","SRPX2"]; axis=["PDPK1","AKT1","MTOR"]
ref=pd.read_csv(J+"M/M3/fig8B_scatter_spearman_long.csv")
fig=new_fig(176)
for i,g in enumerate(genes):
    for j,a_ in enumerate(axis):
        ax=ax_mm(fig,16+j*55,8+i*28,40,20); x=d[a_].values; y=d[g].values
        ax.scatter(x,y,s=1.5,color="#444444",alpha=0.55,linewidths=0,rasterized=True); m,b=np.polyfit(x,y,1); xs=np.array([x.min(),x.max()]); ax.plot(xs,m*xs+b,color="#C0392B",lw=0.8)
        rho,p=spearmanr(x,y); r0=ref[(ref.core_gene==g)&(ref.axis_gene==a_)].rho.values[0]; assert abs(rho-r0)<1e-6,(g,a_,rho,r0)
        ax.set_title(f"rho = {rho:.3f}, p = {p:.1e}",fontsize=8,pad=1.5,loc="left")
        if j==0: ax.set_ylabel(g,fontsize=8,labelpad=1)
        if i==5: ax.set_xlabel(a_,fontsize=8,labelpad=1)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
letter(fig,1,1,"a")
out=V+"FigS9"; save(fig,out); print(check(out+".pdf",V+"FigS9_check.json"))
