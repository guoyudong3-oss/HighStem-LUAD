# Fig2 A-F from the cell-level table and malignant-cell UMAP coordinates
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
from scipy.stats import gaussian_kde, pearsonr
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/Fig2/"; D=V+"data/"
cl=pd.read_csv(D+"fig2_cell_level_data.csv"); um=pd.read_csv(D+"umap_coords.csv"); df=cl.merge(um,on="cell_id")
assert len(df)==13425, len(df)
df["grp"]=pd.Categorical(df.stem_group,["LowStem","DTStem","HighStem"],ordered=True)
cols={"LowStem":"#56B1A6","DTStem":"#C7C7C7","HighStem":"#C92D45"}
keep=df[df.keep_for_umap_plot.astype(str).str.upper()=="TRUE"]
cmap_sp=plt.get_cmap("Spectral_r"); rng=np.random.default_rng(3)
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
def corner_axes(a,xl="UMAP_1",yl="UMAP_2"):
    a.set_xticks([]); a.set_yticks([])
    for s in a.spines.values(): s.set_visible(False)
    x0,x1=a.get_xlim(); y0,y1=a.get_ylim(); L=0.38*(x1-x0)
    a.annotate("",xy=(x0+L,y0),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>",lw=0.7,color="black",shrinkA=0,shrinkB=0,mutation_scale=6))
    a.annotate("",xy=(x0,y0+L),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>",lw=0.7,color="black",shrinkA=0,shrinkB=0,mutation_scale=6))
    a.text(x0+0.66*L,y0-0.03*(y1-y0),xl,fontsize=8,ha="center",va="top"); a.text(x0-0.02*(x1-x0),y0+L/2,yl,fontsize=8,ha="right",va="center",rotation=90)
fig=new_fig(136)
X=[12,72,132]; Y=[8,72]; AW=40; AH=46
# A: CytoTRACE UMAP
a=ax_mm(fig,X[0],Y[0],AW,AH); o=keep.sample(frac=1,random_state=3)
sc=a.scatter(o.UMAP_1,o.UMAP_2,c=o.cytotrace_score,cmap=cmap_sp,vmin=0,vmax=1,s=0.5,linewidths=0,rasterized=True)
a.set_aspect("equal"); a.set_xticks([-5,0,5]); a.set_yticks([-4,0,4]); a.set_xlabel("Component1",fontsize=8,labelpad=1); a.set_ylabel("Component2",fontsize=8,labelpad=1); clean(a)
cax=ax_mm(fig,X[0]+AW+2,Y[0]+10,2.4,26); cb=fig.colorbar(sc,cax=cax,ticks=[0,0.5,1]); cb.ax.tick_params(labelsize=8,length=2,width=0.5); cb.outline.set_linewidth(0.4)
fig.text((X[0]+AW+1.5)/W_MM,1-(Y[0]+6.5)/fig._h_mm,"Less diff.",fontsize=8,ha="left",va="center"); fig.text((X[0]+AW+1.5)/W_MM,1-(Y[0]+40)/fig._h_mm,"More diff.",fontsize=8,ha="left",va="center")
fig.text((X[0]+AH/2)/W_MM,1-(Y[0]-2.2)/fig._h_mm,"CytoTRACE",fontsize=8,fontweight="bold",ha="center",va="center")
# B: rank-percentile histogram
b=ax_mm(fig,X[1],Y[0],AW+8,AH); pct=df.cytotrace_percentile.values
h,e=np.histogram(pct,bins=85,range=(pct.min(),pct.max())); w=np.diff(e); dens=h/(h.sum()*w)
b.bar(e[:-1],dens,width=w,align="edge",color="#EBEBEB",edgecolor="#BFBFBF",linewidth=0.2)
kd=gaussian_kde(pct); xs=np.linspace(0,1,300); b.plot(xs,kd(xs),color="#CC58D8",ls="--",lw=1.0)
q1d=0.25; q3d=0.75
for q in (q1d,q3d): b.plot([q,q],[0,1.18],color="#C0392B",ls="--",lw=0.9)
n={k:int((df.grp==k).sum()) for k in cols}
b.text(0.125,1.5,f"LowStem\n{n['LowStem']:,} cells",fontsize=8,ha="center",va="top",linespacing=0.95); b.text(0.5,1.5,f"DTStem\n{n['DTStem']:,} cells",fontsize=8,ha="center",va="top",linespacing=0.95); b.text(0.875,1.5,f"HighStem\n{n['HighStem']:,} cells",fontsize=8,ha="center",va="top",linespacing=0.95)
b.set_xlim(-0.03,1.03); b.set_ylim(0,1.55); b.set_xticks([0,0.25,0.5,0.75,1]); b.set_xticklabels(["0","0.25","0.5","0.75","1"]); b.set_yticks([0,0.5,1.0]); b.set_xlabel("Rank percentile of CytoTRACE score",fontsize=8,labelpad=1); b.set_ylabel("Density",fontsize=8,labelpad=1); clean(b)
b.set_title("Stem score distribution",fontsize=8,fontweight="bold",pad=2)
# C: stemness groups on UMAP
c=ax_mm(fig,X[2],Y[0],AW,AH); o=keep.sample(frac=1,random_state=4)
c.scatter(o.UMAP_1,o.UMAP_2,c=[cols[k] for k in o.grp],s=0.5,linewidths=0,rasterized=True); c.set_aspect("equal"); c.set_xlim(c.get_xlim()); c.set_ylim(c.get_ylim()); corner_axes(c)
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=4.5,color=cols[k],label=k) for k in cols],loc="upper left",ncol=3,bbox_to_anchor=((X[2]-10)/W_MM,1-(Y[0]+AH+5)/fig._h_mm),frameon=False,fontsize=8,handletextpad=-0.3,columnspacing=0.4,borderaxespad=0)
fig.text((X[2]+AH/2)/W_MM,1-(Y[0]-2.2)/fig._h_mm,"Malignant cells",fontsize=8,fontweight="bold",ha="center",va="center")
# D, E: raincloud plots
def rain(a,col,ylab,ylim,yt,trim,adj):
    for i,g in enumerate(cols,1):
        v=df[df.grp==g][col].values
        kde=gaussian_kde(v,bw_method=adj); ys=np.linspace(v.min() if trim else v.min()-0.1*np.ptp(v),v.max() if trim else v.max()+0.1*np.ptp(v),200); d=kde(ys); d=d/d.max()*0.72*0.5
        a.fill_betweenx(ys,i+0.04,i+0.04+d,color=cols[g],alpha=0.62,lw=0)
        bp=a.boxplot([v],positions=[i],widths=0.14,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=0.8),boxprops=dict(lw=0.6,facecolor=cols[g]),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6))
        a.scatter(np.full(len(v),0.56)+rng.uniform(-0.025,0.025,len(v)),v,s=0.8,color=cols[g],linewidths=0,rasterized=True)
    ymax=df[col].max(); st=(ylim[1]-ylim[0])*0.065; tk=st*0.35
    for (x0,x1,yy) in [(1.02,1.98,ymax+st),(2.02,2.98,ymax+2*st),(1.0,3.0,ymax+3.2*st)]:
        a.plot([x0,x0,x1,x1],[yy-tk,yy,yy,yy-tk],color="black",lw=0.6,clip_on=False); a.text((x0+x1)/2,yy-tk*0.4,"***",fontsize=8,ha="center",va="bottom")
    a.set_xlim(0.4,3.65); a.set_ylim(*ylim); a.set_yticks(yt); a.set_xticks([1,2,3]); a.set_xticklabels(list(cols),rotation=45,ha="right",rotation_mode="anchor"); a.set_ylabel(ylab,fontsize=8,labelpad=1); clean(a)
d_=ax_mm(fig,X[0],Y[1],AW+4,AH); rain(d_,"cytotrace_score","CytoTRACE score",(0,1.22),[0,0.25,0.5,0.75,1.0],False,0.75)
e_=ax_mm(fig,X[1]+4,Y[1],AW+4,AH); rain(e_,"trs_score","TRS score",(0.15,1.32),[0.5,1.0],True,None)
# F: scatter with marginal histograms
f=ax_mm(fig,X[2]+2,Y[1]+8,AW-6,AH-8); r,p=pearsonr(df.cytotrace_score,df.trs_score); o=df.sample(frac=1,random_state=5)
f.scatter(o.cytotrace_score,o.trs_score,s=1.2,color="#555555",linewidths=0,rasterized=True)
m,bb=np.polyfit(df.cytotrace_score,df.trs_score,1); xx=np.array([0,1]); f.plot(xx,m*xx+bb,color="blue",lw=1.0)
f.set_xlim(0,1.0); f.set_ylim(0.1,1.05); f.set_xlabel("CytoTRACE score",fontsize=8,labelpad=1); f.set_ylabel("TRS score",fontsize=8,labelpad=1); clean(f); f.set_xticks([0,0.5,1.0]); f.set_yticks([0.2,0.6,1.0])
ft=ax_mm(fig,X[2]+2,Y[1],AW-6,7.5); ft.hist(df.cytotrace_score,bins=32,color="#16A085",edgecolor="#0E6655",linewidth=0.2); ft.set_xlim(0,1.0); ft.axis("off")
fr=ax_mm(fig,X[2]+2+AW-6+0.8,Y[1]+8,7.5,AH-8); fr.hist(df.trs_score,bins=32,orientation="horizontal",color="#E67E22",edgecolor="#A04000",linewidth=0.2); fr.set_ylim(0.1,1.05); fr.axis("off")
fig.text((X[2]+2)/W_MM,1-(Y[1]-3)/fig._h_mm,f"Pearson r = {r:.2f}, 95% CI 0.60–0.62, n = 13,425, p < 2.2e-16",fontsize=8,ha="left",va="center") if False else None
zz=np.arctanh(r); se=1/np.sqrt(len(df)-3); lo,hi=np.tanh(zz-1.959964*se),np.tanh(zz+1.959964*se); fig.text((X[2]+2)/W_MM,1-(Y[1]-5.5)/fig._h_mm,f"r = {r:.3f} (95% CI {lo:.3f}–{hi:.3f})",fontsize=8,va="center"); fig.text((X[2]+2)/W_MM,1-(Y[1]-1.5)/fig._h_mm,"n = 13,425, p < 2.2e-16",fontsize=8,va="center")
for t,x,y in [("a",1,1.5),("b",X[1]-9.5,1.5),("c",X[2]-9.5,1.5),("d",1,Y[1]-6.5),("e",X[1]-5.5,Y[1]-6.5),("f",X[2]-9.5,Y[1]-6.5)]: letter(fig,x,y,t)
out=V+"Fig2_main_v3"; save(fig,out); print(check(out+".pdf",V+"Fig2_main_check.json"))
