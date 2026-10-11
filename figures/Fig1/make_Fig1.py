# Fig1: A UMAP clusters, B cell types, C marker dot plot (genes on the x axis)
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
import matplotlib.patheffects as pe
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/Fig1/"; D=V+"data/"
cells=pd.read_csv(D+"fig1_cells_umap_clusters_celltypes.csv",dtype={"seurat_clusters":str}); dot=pd.read_csv(D+"fig1C_dot_data.csv"); genes=[g.strip() for g in open(D+"fig1C_marker_genes.txt")]
npg=["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85"]
def ramp(cols,n):
    rgb=np.array([mpl.colors.to_rgb(c) for c in cols]); x=np.linspace(0,1,len(cols)); xi=np.linspace(0,1,n); return [tuple(np.interp(xi,x,rgb[:,k])[i] for k in range(3)) for i in range(n)]
cl=sorted(cells.seurat_clusters.unique()); ccol=dict(zip(cl,ramp(npg,len(cl))))
pal={"NK":"#4482B7","T":"#F17E24","Plasma":"#4AB035","Epi":"#DA2F2D","Fib":"#9367C6","B":"#8A554B","Endo":"#DE7ACE","Mac":"#B7C838","Mon":"#B6CC63","Mural":"#62D4D2","Cyc":"#B1CBE7","Mast":"#EEB778","MDC":"#A3F092","Glial":"#F39396","Other":"#CCCCCC"}
rng=np.random.default_rng(1); order=rng.permutation(len(cells)); c=cells.iloc[order]
fig=new_fig(172)
# A
axA=ax_mm(fig,9,4,78,78); axA.scatter(c.UMAP_1,c.UMAP_2,s=0.15,c=[ccol[k] for k in c.seurat_clusters],linewidths=0,rasterized=True)
med=cells.groupby("seurat_clusters")[["UMAP_1","UMAP_2"]].median()
for k,(x,y) in med.iterrows(): axA.text(x,y,k,fontsize=8,ha="center",va="center",fontweight="bold",path_effects=[pe.withStroke(linewidth=1.6,foreground="white")])
# B
axB=ax_mm(fig,97,4,62,78); axB.scatter(c.UMAP_1,c.UMAP_2,s=0.15,c=[pal[k] for k in c.celltype_display],linewidths=0,rasterized=True)
for a in (axA,axB):
    a.set_xticks([]); a.set_yticks([]); a.set_xlabel("UMAP_1",fontsize=8,labelpad=1); a.set_ylabel("UMAP_2",fontsize=8,labelpad=1)
    for s in a.spines.values(): s.set_linewidth(0.6)
lev=["NK","T","B","Plasma","Mon","Mac","Mast","Fib","Mural","Endo","Epi","Cyc","Glial","MDC","Other"]
from matplotlib.lines import Line2D
h=[Line2D([],[],marker="o",ls="",markersize=4.5,color=pal[k],label=k) for k in lev]
fig.legend(handles=h,loc="center left",bbox_to_anchor=(161/W_MM,1-43/fig._h_mm),frameon=False,fontsize=8,handletextpad=0.1,labelspacing=0.35)
letter(fig,1,1.5,"a"); letter(fig,91,1.5,"b")
# C: dot plot (x = genes, y = cell types)
ctypes=["NK","T","Plasma","Epi","Fib","B","Endo","Mac","Mon","Mural","Cyc","Mast","MDC","Glial"]
q=np.nanquantile(dot["avg.exp.scaled"],[0.02,0.98]); v=np.clip((dot["avg.exp.scaled"]-q[0])/(q[1]-q[0]),0,1)**0.7; dot["v"]=v; dot["pct"]=dot["pct.exp"].clip(0,75)
xi={g:i for i,g in enumerate(genes)}; yi={t:i for i,t in enumerate(ctypes[::-1])}
axC=ax_mm(fig,19,92,141,55)
cmap=mpl.colors.LinearSegmentedColormap.from_list("f1",["#2166AC","#F7F7F7","#67001F"])
sz=lambda p:(p/75*1.0)**1.0*9.0
sc=axC.scatter([xi[g] for g in dot["features.plot"]],[yi[t] for t in dot["id"]],s=[max(sz(p),0.05) for p in dot.pct],c=dot.v,cmap=cmap,vmin=0,vmax=1,linewidths=0)
axC.set_xlim(-0.8,len(genes)-0.2); axC.set_ylim(-0.7,len(ctypes)-0.3)
axC.set_xticks(range(len(genes))); axC.set_xticklabels(genes,rotation=90,fontsize=8); axC.set_yticks(range(len(ctypes))); axC.set_yticklabels(ctypes[::-1],fontsize=8)
axC.tick_params(length=0,pad=1.5)
for s in axC.spines.values(): s.set_linewidth(0.6)
# colour strip marking the marker-gene group of each cell type
grp={"NK":["NKG7","KLRD1","GNLY"],"T":["CD8A","CD3D","CD3E","CCL5"],"Plasma":["JCHAIN","IGKC","IGHG1"],"Epi":["EPCAM","CDH1","KRT7","PGC","PGA3","MUC5AC"],"Fib":["PDGFRA","LUM","DCN"],"B":["CD79A","CD79B","MS4A1"],"Endo":["VWF","PECAM1"],"Mac":["C1QA","MARCO","C1QB","CD68"],"Mon":["LYZ","FCN1","S100A9","S100A8"],"Mural":["RGS5","PDGFRB","ACTA2"],"Cyc":["MKI67","TOP2A","CDK1"],"Mast":["CPA3","GATA2","KIT"],"MDC":["CD1C","FCER1A","CLEC9A"],"Glial":["CDH19","PLP1","SOX10"]}
for t,gs in grp.items():
    x0=min(xi[g] for g in gs if g in xi)-0.45; x1=max(xi[g] for g in gs if g in xi)+0.45
    axC.add_patch(plt.Rectangle((x0,len(ctypes)-0.3+0.05),x1-x0,0.35,color=pal[t],clip_on=False,lw=0))
axC.set_ylim(-0.7,len(ctypes)+0.15)
# legends below the dot plot
cax=ax_mm(fig,60,165,32,2.6); cb=fig.colorbar(sc,cax=cax,orientation="horizontal",ticks=[0,0.5,1]); cb.ax.tick_params(labelsize=8,length=2,width=0.5,pad=1); cb.outline.set_linewidth(0.4)
fig.text(19/W_MM,1-166.3/fig._h_mm,"Relative expression",fontsize=8,va="center",ha="left")
fig.text(105/W_MM,1-166.3/fig._h_mm,"Percent expressed",fontsize=8,va="center",ha="left")
for k,p in enumerate([25,50,75]):
    x=140+k*13; fig.add_artist(mpl.lines.Line2D([x/W_MM],[1-166.3/fig._h_mm],marker="o",ms=np.sqrt(sz(p)),color="black",ls="",transform=fig.transFigure)); fig.text((x+2.2)/W_MM,1-166.3/fig._h_mm,str(p),fontsize=8,va="center")
letter(fig,1,87,"c")
out=V+"Fig1_main"; save(fig,out); print(check(out+".pdf",V+"Fig1_main_check.json"))
