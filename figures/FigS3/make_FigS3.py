# Fig S3: A circle plots and B role heatmaps (embedded R output); C, D ligand-receptor bubbles and E metabolic pathways from tables
import sys, textwrap; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/FigS3/"; t=J+"J4/main_1200/tables/"
short={"B lymphocytes":"B","DTStem":"DTStem","Endothelial cells":"Endo","Fibroblasts":"Fib","HighStem":"HighStem","LowStem":"LowStem","MAST cells":"Mast","Myeloid cells":"Myeloid","NK cells":"NK","NormalEpi":"NormalEpi","Oligodendrocytes":"Oligodend.","T lymphocytes":"T","T/NK cells":"T/NK","Undetermined":"Undeterm."}
NLR=20
def lr_table(f,partner_col):
    d=pd.read_csv(t+f); d=d[d.pval<=0.01].copy()   # as in 18b: p <= 0.01
    keep=d.groupby("interaction_name_2").prob.max().sort_values(ascending=False).head(NLR).index.tolist()   # as in 18b: ranked by maximum communication probability
    d=d[d.interaction_name_2.isin(keep)].copy(); d["partner"]=d[partner_col].map(short); return d,keep
dC,kC=lr_table("fig3F_source_table_highstem_outgoing.csv","target"); dD,kD=lr_table("fig3G_source_table_highstem_incoming.csv","source")
assert (dC.source=="HighStem").all() and (dD.target=="HighStem").all()
vmax=max(dC.prob.max(),dD.prob.max())*1e3
fig=new_fig(230)
partners=["B","DTStem","Endo","Fib","LowStem","Mast","Myeloid","NK","NormalEpi","Oligodend.","T","T/NK","Undeterm."]
cmap=mpl.colors.LinearSegmentedColormap.from_list("cp",["#2C7BB6","#66C2A5","#FEE08B","#F46D43","#D73027"])
def bubble(x,y,w,h,d,keep,title):
    ax=ax_mm(fig,x,y,w,h); xi={p:i for i,p in enumerate(partners)}; yi={l:len(keep)-1-i for i,l in enumerate(keep)}
    sc=ax.scatter([xi[p] for p in d.partner],[yi[l] for l in d.interaction_name_2],c=d.prob*1e3,cmap=cmap,vmin=0,vmax=vmax,s=11,linewidths=0)
    ax.set_xticks(range(len(partners))); ax.set_xticklabels(partners,rotation=90); ax.set_yticks(range(len(keep))); ax.set_yticklabels(keep[::-1])
    ax.set_xlim(-0.6,len(partners)-0.4); ax.set_ylim(-0.6,len(keep)-0.4); ax.tick_params(length=1.5,pad=1); ax.grid(True,color="#ECECEC",lw=0.3); ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_linewidth(0.5); s.set_color("#777777")
    ax.set_title(title,fontsize=8,fontweight="bold",pad=2); return sc
Y3=154
sc=bubble(38,Y3,38,57,dC,kC,"HighStem → partner (outgoing)")
bubble(126,Y3,38,57,dD,kD,"Partner → HighStem (incoming)")
cax=ax_mm(fig,167,Y3+6,2.4,26); cb=fig.colorbar(sc,cax=cax); cb.ax.tick_params(labelsize=8,length=2,width=0.5); cb.outline.set_linewidth(0.4)
fig.text(165/W_MM,1-(Y3+38)/fig._h_mm,"Prob.\n(×0.001)",fontsize=8,va="center",linespacing=0.95)
fig.text(163/W_MM,1-(Y3+46)/fig._h_mm,"all p < 0.01",fontsize=8,va="center")
# E: metabolic pathways ranked by HighStem advantage (as in 18b); top 4 and bottom 4 shown
mi=pd.read_csv(t+"fig3I_metabolic_bubble_values.csv"); w_=mi.pivot_table(index="pathway",columns="stem_group",values="score"); w_["adv"]=w_.HighStem-w_[["DTStem","LowStem"]].mean(1); w_=w_.sort_values("adv",ascending=False)
sel=list(w_.index[:4])+list(w_.index[-4:]); grp=["LowStem","DTStem","HighStem"]
e=ax_mm(fig,146,5,13,46); vals=mi.set_index(["pathway","stem_group"]).value
cm2=mpl.colors.LinearSegmentedColormap.from_list("m",["#67C5CF","#F0E442","#E53935"])
for r,pw in enumerate(sel[::-1]):
    for c,g in enumerate(grp):
        v=vals[(pw,g)]; e.scatter(c,r,s=4+v*28,color=cm2(v),linewidths=0)
e.axhline(3.5,color="#999999",lw=0.5,ls="--")
e.set_xticks(range(3)); e.set_xticklabels(grp,rotation=90); e.set_yticks(range(len(sel))); e.set_yticklabels([textwrap.fill(s,30) for s in sel[::-1]],linespacing=0.9); e.set_xlim(-0.6,2.6); e.set_ylim(-0.6,len(sel)-0.4)
e.tick_params(length=1.5,pad=1); [s.set_linewidth(0.5) for s in e.spines.values()]
cax2=ax_mm(fig,166,12,2.4,20); cb2=fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0,1),cmap=cm2),cax=cax2,ticks=[0,0.5,1]); cb2.ax.tick_params(labelsize=8,length=2,width=0.5); cb2.outline.set_linewidth(0.4)
fig.text(164/W_MM,1-8/fig._h_mm,"Value",fontsize=8,va="center")
for tt,x,y in [("a",1,1),("b",1,63),("c",1,150),("d",89,150),("e",104,1)]: letter(fig,x,y,tt)
tmp=V+"FigS3_v2_base.pdf"; fig.savefig(tmp,transparent=True)
embed_pdf(tmp,V+"FigS3_v2.pdf",[(V+"panels/S3A_circles_v2.pdf",0,2),(V+"panels/S3B_role_heatmaps_v2.pdf",0,64)])
pd.DataFrame({"panel":["C"]*NLR+["D"]*NLR,"rank":list(range(1,NLR+1))*2,"LR":kC+kD}).to_csv(V+"FigS3_v2_CD_LR_shown.csv",index=False)
pd.DataFrame({"pathway":sel,"HighStem_advantage":w_.loc[sel,"adv"].values}).to_csv(V+"FigS3_v2_E_pathways_shown.csv",index=False)
print(check(V+"FigS3_v2.pdf",V+"FigS3_v2_check.json"))
