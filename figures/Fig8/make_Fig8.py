# Fig8: A pathway scores by stemness group, B gene x PI3K-axis Spearman rho, C six-gene score vs PI3K-axis genes
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/Fig8/"
z=pd.read_csv(J+"J4/main_1200/tables/fig3H_pathway_activity_group_zscore.csv",index_col=0)
foc=["PI3K","Hypoxia","NFkB","JAK-STAT"]; grp=["LowStem","DTStem","HighStem"]; Z=z.loc[foc,grp]
hm=pd.read_csv(J+"L/L2/cache/fig8_core_vs_pi3k_axis_heatmap_spearman.csv"); sg=pd.read_csv(J+"L/L2/cache/fig8_signature_vs_pi3k_axis_spearman.csv")
genes=["TMSB10","S100A16","HLA-C","S100A14","CKLF","SRPX2"]; axis=["PIK3CA","PIK3CB","PDPK1","AKT1","MTOR"]
hm["core_gene"]=hm.core_gene.str.replace(".","-",regex=False); R=hm.pivot(index="core_gene",columns="axis_gene",values="rho").loc[genes,axis]; P=hm.pivot(index="core_gene",columns="axis_gene",values="p_value").loc[genes,axis]
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
fig=new_fig(130)
# A
a=ax_mm(fig,26,12,50,28); cm=mpl.colors.LinearSegmentedColormap.from_list("t",["#38807F","#F0F0F0","#AB2329"]); zm=np.abs(Z.values).max()
im=a.imshow(Z.values,cmap=cm,vmin=-zm,vmax=zm,aspect="auto")
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]): a.text(j,i,f"{Z.values[i,j]:.2f}",ha="center",va="center",fontsize=8,color="white" if abs(Z.values[i,j])>0.7*zm else "black")
a.set_xticks(range(3)); a.set_xticklabels(grp,rotation=40,ha="right",rotation_mode="anchor"); a.set_yticks(range(4)); a.set_yticklabels(foc); a.tick_params(length=0)
for s in a.spines.values(): s.set_visible(False)
cax=ax_mm(fig,80,14,2.4,24); cb=fig.colorbar(im,cax=cax); cb.ax.tick_params(labelsize=8,length=2,width=0.5); cb.outline.set_linewidth(0.4); cb.set_label("Standardized score",fontsize=8,labelpad=2)
fig.text(26/W_MM,1-8/fig._h_mm,"Pathway activity in malignant cells (13,425)",fontsize=8,fontweight="bold",va="center") if False else None
# B: heatmap
b=ax_mm(fig,112,12,50,40); bm=mpl.colors.LinearSegmentedColormap.from_list("rho",["#FFFFFF","#FDD0A2","#E6550D","#7F2704"]); im2=b.imshow(R.values,cmap=bm,vmin=0.1,vmax=0.55,aspect="auto")
for i in range(6):
    for j in range(5): b.text(j,i,f"{R.values[i,j]:.2f}",ha="center",va="center",fontsize=8,color="white" if R.values[i,j]>0.42 else "black")
b.set_xticks(range(5)); b.set_xticklabels(axis,rotation=40,ha="right",rotation_mode="anchor"); b.set_yticks(range(6)); b.set_yticklabels(genes); b.tick_params(length=0)
for s in b.spines.values(): s.set_visible(False)
cax=ax_mm(fig,164,14,2.4,24); cb=fig.colorbar(im2,cax=cax,ticks=[0.1,0.3,0.5]); cb.ax.tick_params(labelsize=8,length=2,width=0.5); cb.outline.set_linewidth(0.4); cb.set_label("Spearman rho",fontsize=8,labelpad=2); (lambda *a,**k:None)(158/W_MM,1-11.5/fig._h_mm,"rho",fontsize=8,va="bottom")
# C: lollipop
c=ax_mm(fig,26,86,70,30); sg=sg.sort_values("rho"); y=np.arange(len(sg))
c.hlines(y,0,sg.rho,color="#BBBBBB",lw=0.8); c.scatter(sg.rho,y,s=26,c=sg.rho,cmap=bm,vmin=0.1,vmax=0.6,edgecolor="black",linewidth=0.3,zorder=3)
for yi,(r,p) in enumerate(zip(sg.rho,sg.p_value)): c.text(r+0.03,yi,f"{r:.3f}  (p = {p:.1e})",fontsize=8,va="center")
c.set_yticks(y); c.set_yticklabels(sg.gene); c.set_xlim(0,1.05); c.set_xticks([0,0.2,0.4,0.6]); c.set_xlabel("Spearman rho (six-gene program score)",fontsize=8,labelpad=1); clean(c); c.spines["left"].set_visible(False); c.tick_params(axis="y",length=0)
for t,x,y_ in [("a",1,1),("b",96,1),("c",1,74)]: letter(fig,x,y_,t)
fig.text(26/W_MM,1-7/fig._h_mm,"Pathway activity, malignant cells",fontsize=8,fontweight="bold",va="center")
fig.text(112/W_MM,1-7/fig._h_mm,"Six genes vs PI3K axis (n = 517)",fontsize=8,fontweight="bold",va="center")
fig.text(26/W_MM,1-80/fig._h_mm,"Six-gene program score vs PI3K-axis genes (n = 517)",fontsize=8,fontweight="bold",va="center")
out=V+"Fig8_main_v2"; save(fig,out); print(check(out+".pdf",V+"Fig8_main_check.json"))
