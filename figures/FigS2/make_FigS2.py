# Fig S2: A CNV heatmap (embedded R output), B CNV score by k cluster, C TRS on UMAP, D stem/TRS overlay with colour key
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
from scipy.stats import gaussian_kde
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/"; D=V+"Fig2/data/"
cl=pd.read_csv(D+"fig2_cell_level_data.csv"); um=pd.read_csv(D+"umap_coords.csv"); df=cl.merge(um,on="cell_id"); keep=df[df.keep_for_umap_plot.astype(str).str.upper()=="TRUE"]
kc=pd.read_csv(D+"k_cluster_all.csv"); kpal={1:"#8DD3C7",2:"#FFFFB3",3:"#B3DE69",4:"#FB8072"}
AH=70.0; fig=new_fig(AH+120)
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
Y0=AH+8
# B: violins
b=ax_mm(fig,16,Y0+4,50,42)
for k in (1,2,3,4):
    v=kc[kc.k_cluster==k].cnv_score_runfinal.values; kd=gaussian_kde(v); ys=np.linspace(v.min()-0.1*np.ptp(v),v.max()+0.1*np.ptp(v),300); d=kd(ys); d=d/d.max()*0.45
    b.fill_betweenx(ys,k-d,k+d,facecolor=kpal[k],edgecolor="black",lw=0.4,alpha=0.85)
    b.boxplot([v],positions=[k],widths=0.2,showfliers=False,patch_artist=True,boxprops=dict(facecolor="white",lw=0.5),medianprops=dict(color="black",lw=0.8),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5))
b.set_xticks([1,2,3,4]); b.set_xlabel("k cluster",fontsize=8,labelpad=1); b.set_ylabel("CNV score",fontsize=8,labelpad=1); clean(b); b.set_ylim(-0.01,0.26)
b.set_title("Kruskal-Wallis, p < 2.2e-16",fontsize=8,pad=2)
# C TRS UMAP
c=ax_mm(fig,84,Y0+4,38,42); o=keep.sample(frac=1,random_state=6)
sc=c.scatter(o.UMAP_1,o.UMAP_2,c=o.trs_score,cmap="viridis",s=0.5,linewidths=0,rasterized=True); c.set_aspect("equal")
def corner(a,xl="UMAP_1",yl="UMAP_2"):
    a.set_xticks([]); a.set_yticks([])
    for s in a.spines.values(): s.set_visible(False)
    x0,x1=a.get_xlim(); y0,y1=a.get_ylim(); L=0.38*(x1-x0)
    for tx,ty in [((x0+L,y0),(x0,y0)),((x0,y0+L),(x0,y0))]: a.annotate("",xy=tx,xytext=ty,arrowprops=dict(arrowstyle="-|>",lw=0.7,color="black",shrinkA=0,shrinkB=0,mutation_scale=6))
    a.text(x0+0.66*L,y0-0.03*(y1-y0),xl,fontsize=8,ha="center",va="top"); a.text(x0-0.02*(x1-x0),y0+L/2,yl,fontsize=8,ha="right",va="center",rotation=90)
c.set_xlim(c.get_xlim()); c.set_ylim(c.get_ylim()); corner(c)
cax=ax_mm(fig,Y0*0+124.5,Y0+12,2.4,26); cb=fig.colorbar(sc,cax=cax,ticks=[0.25,0.5,0.75,1.0]); cb.ax.tick_params(labelsize=8,length=2,width=0.5); cb.outline.set_linewidth(0.4)
fig.text(131/W_MM,1-(Y0+9)/fig._h_mm,"TRS score",fontsize=8,ha="left",va="center")
fig.text(103/W_MM,1-(Y0+1.5)/fig._h_mm,"scPagwas TRS score",fontsize=8,fontweight="bold",ha="center",va="center")
# D: overlays and colour key
Y1=Y0+58
rs=lambda x:(x-x.min())/(x.max()-x.min())
sn=rs(keep.cytotrace_score.values); tn=rs(keep.trs_score.values); grey=np.array([0.8,0.8,0.8])
def comp(col,alpha): a=np.asarray(alpha)[:,None]; return a*np.array(col)[None,:]+(1-a)*grey[None,:]
def amap(x,lo,hi,mx): return mx*(x-x.min())/(x.max()-x.min())
a1=amap(np.maximum(0.08,sn),0,1,0.95); a2=amap(np.maximum(0.08,tn),0,1,0.95); mix=np.stack([sn,tn,np.zeros_like(sn)],1); a3=amap(np.maximum(0.10,np.maximum(sn,tn)),0,1,0.92)
cols1=comp(mpl.colors.to_rgb("#EE0000"),a1); cols2=comp(mpl.colors.to_rgb("#6EDB78"),a2); cols3=a3[:,None]*mix+(1-a3[:,None])*grey[None,:]
order=np.random.default_rng(7).permutation(len(keep))
PW=36; xs=[10,52,94]
for x,cc,tt in zip(xs,[cols1,cols2,cols3],["Stem score","TRS score","Stem + TRS"]):
    a=ax_mm(fig,x,Y1+6,PW,PW); a.scatter(keep.UMAP_1.values[order],keep.UMAP_2.values[order],c=cc[order],s=0.5,linewidths=0,rasterized=True); a.set_aspect("equal"); a.set_xlim(a.get_xlim()); a.set_ylim(a.get_ylim()); corner(a)
    a.set_title(tt,fontsize=8,fontweight="bold",pad=2)
ks=ax_mm(fig,138,Y1+6,32,PW)
G=np.zeros((10,10,3)); 
for i in range(10):
    for j in range(10): G[j,i]=[(i+1)/10,(j+1)/10,0]
ks.imshow(G,origin="lower",extent=(0.5,10.5,0.5,10.5),aspect="auto"); ks.set_xticks([2,4,6,8,10]); ks.set_yticks([2.5,5,7.5,10]); ks.set_yticklabels(["2.5","5.0","7.5","10"]); ks.set_xlabel("Stem score (1–10)",fontsize=8,labelpad=1); ks.set_ylabel("TRS score (1–10)",fontsize=8,labelpad=1)
ks.set_title("Colour key",fontsize=8,fontweight="bold",pad=2); clean(ks)
for t,x,y in [("b",1,Y0-2),("c",74,Y0-2),("d",1,Y1-2)]: letter(fig,x,y,t)
letter(fig,1,1.5,"a")
tmp=V+"FigS2/FigS2_v2_lower.pdf"; fig.savefig(tmp,transparent=True)
embed_pdf(tmp,V+"FigS2/FigS2_v2.pdf",[(V+"FigS2/A/out/FigS2A_cnv_heatmap_v2.pdf",0,0)])
print(check(V+"FigS2/FigS2_v2.pdf",V+"FigS2/FigS2_v2_check.json"))
