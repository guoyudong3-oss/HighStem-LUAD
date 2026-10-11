# Fig3 A-D: outgoing/incoming strength, partner-level pair counts, pathway heatmap (embedded R output)
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
from scipy.stats import wilcoxon
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/Fig3/"; t=J+"J4/main_1200/tables/"
cm=pd.read_csv(t+"fig3_cellchat_count_matrix.csv",index_col=0); w=pd.read_csv(t+"fig3_cellchat_weight_matrix.csv",index_col=0)
gcol={"B lymphocytes":"#E41A1C","DTStem":"#377EB8","Endothelial cells":"#4DAF4A","Fibroblasts":"#984EA3","HighStem":"#F18F01","LowStem":"#F781BF","MAST cells":"#B39BC8","Myeloid cells":"#A65628","NK cells":"#56B9E9","NormalEpi":"#1F2A7C","T lymphocytes":"#1B9E77","T/NK cells":"#9BCF6E","Undetermined":"#D9B300","Oligodendrocytes":"#DDA08A"}
fig=new_fig(184)
# A: scatter
W2=w.copy().values; np.fill_diagonal(W2,0); C2=cm.copy().values; np.fill_diagonal(C2,0)
d=pd.DataFrame({"label":w.index,"x":W2.sum(1),"y":W2.sum(0),"count":C2.sum(1)+C2.sum(0)}); d=d[(d.x>0)&(d.y>0)].reset_index(drop=True)
qx,qy=d.x.quantile(0.95),d.y.quantile(0.95); d["xp"]=d.x.clip(upper=qx); d["yp"]=d.y.clip(upper=qy)
a=ax_mm(fig,14,9,70,66)
a.scatter(d.xp,d.yp,s=(d["count"]/d["count"].max())**0.5*70+10,c=[gcol[l] for l in d.label],linewidths=0,zorder=3)
# label positions (data coordinates; alignment; leader line)
lab={"HighStem":(0.127,0.103,"left",0),"Myeloid cells":(0.065,0.1095,"center",0),"Fibroblasts":(0.154,0.0545,"center",0),"DTStem":(0.094,0.073,"left",0),"NormalEpi":(0.05,0.0795,"center",0),
 "T/NK cells":(0.0,0.0795,"right",1),"T lymphocytes":(0.0,0.0715,"right",1),"NK cells":(0.0,0.0640,"right",1),"MAST cells":(0.0,0.0570,"right",1),"B lymphocytes":(0.0,0.0500,"right",1),
 "Undetermined":(0.0,0.0425,"right",1),"Endothelial cells":(0.075,0.0415,"center",0),"LowStem":(0.055,0.0285,"center",0),"Oligodendrocytes":(0.052,0.0075,"center",0)}
for _,r in d.iterrows():
    x,y,ha,ln=lab[r.label]; a.text(x,y,r.label,fontsize=8,ha=ha,va="center",zorder=4)
    if ln: a.plot([x+0.002,r.xp],[y,r.yp],color="#999999",lw=0.4,zorder=2)
a.set_xlim(-0.10,0.19); a.set_ylim(0.0,0.125); a.set_xticks([0.05,0.10,0.15]); a.set_yticks([0.025,0.05,0.075,0.10]); a.set_xlabel("Outgoing interaction strength",fontsize=8,labelpad=1); a.set_ylabel("Incoming interaction strength",fontsize=8,labelpad=1)
a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
for k,v in enumerate([200,400,600]):
    fig.add_artist(mpl.lines.Line2D([(66+3)/W_MM],[1-(63+k*5.2)/fig._h_mm],marker="o",ms=np.sqrt((v/d["count"].max())**0.5*70+10),color="black",ls="",transform=fig.transFigure)); fig.text(73/W_MM,1-(63+k*5.2)/fig._h_mm,str(v),fontsize=8,va="center")
fig.text(66/W_MM,1-57.5/fig._h_mm,"Count",fontsize=8,va="center")
# B, C: partner-level bars (HighStem vs LowStem / DTStem)
D0=pd.read_csv(t+"fig3D_high_vs_low_interaction_pairs.csv"); partners=D0.partner.tolist()
show={"NormalEpi":"NormalEpi","Myeloid cells":"Myeloid","Fibroblasts":"Fib","Endothelial cells":"Endo","Undetermined":"Undeterm.","NK cells":"NK","Oligodendrocytes":"Oligodend.","T lymphocytes":"T","T/NK cells":"T/NK","B lymphocytes":"B","MAST cells":"Mast"}
for p in partners: assert p in show,p
def df_for(ref):
    x=pd.DataFrame({"partner":partners,"high":cm.loc["HighStem",partners].values,"ref":cm.loc[ref,partners].values}); x["log2"]=np.log2((x.high+1)/(x.ref+1)); return x.sort_values("high",ascending=False).reset_index(drop=True)
dD,dE=df_for("LowStem"),df_for("DTStem")
def draw(x0,y0,dd,ref):
    ax=ax_mm(fig,x0,y0,60,56); n=len(dd); y=np.arange(n)[::-1]; bmax=dd.high.max(); xl=-dd.ref.max()*1.3; d0=bmax*1.10; d1=bmax*1.50; xr=d1*1.04; lmax=max(2.0,np.ceil(dd.log2.max()*2)/2)
    ax.barh(y,-dd.ref,color="#3A8080",height=0.78); ax.barh(y,dd.high,color="#AB2328",height=0.78)
    for yi,r,h in zip(y,dd.ref,dd.high):
        ax.text(-r+0.03*bmax,yi,str(int(r)),color="white",va="center",ha="left",fontsize=8,fontweight="bold"); ax.text(h-0.02*bmax,yi,str(int(h)),color="white",va="center",ha="right",fontsize=8,fontweight="bold")
    xd=d0+(np.clip(dd.log2,0,lmax)/lmax)*(d1-d0); ax.axvspan(d0-0.02*bmax,xr,color="#F1F1F1",zorder=0)
    ax.scatter(xd,y,s=14,color=np.where(dd.high>dd.ref,"#CC0000","black"),zorder=3,linewidths=0)
    ax.axvline(0,color="grey",lw=0.6); ax.set_xlim(xl,xr); ax.set_ylim(-0.7,n-0.3); ax.set_yticks(y); ax.set_yticklabels([show[p] for p in dd.partner],fontsize=8); ax.yaxis.tick_right(); ax.tick_params(axis="y",length=0,pad=2)
    step=20 if bmax>60 else 10; ax.set_xticks(np.arange(0,bmax*1.05,step)); ax.set_xlabel("# of interaction pairs",fontsize=8,x=0.36,labelpad=1)
    top=ax.secondary_xaxis("top"); tk=np.arange(0,lmax+1e-9,1.0); top.set_xticks(d0+tk/lmax*(d1-d0)); top.set_xticklabels([f"{v:g}" for v in tk],fontsize=8); top.set_xlabel(f"log2 ratio (HighStem / {ref})",fontsize=8,x=0.74,labelpad=1)
    for sp in ("top","right","left","bottom"): ax.spines[sp].set_linewidth(0.5)
    return ax
draw(98,13,dD,"LowStem"); draw(18,96,dE,"DTStem")
# D: heatmap (R output)
letter(fig,1,1.5,"a"); letter(fig,92,1.5,"b"); letter(fig,1,88,"c"); letter(fig,92,88,"d")
tmp=V+"Fig3_base.pdf"; fig.savefig(tmp,transparent=True)
embed_pdf(tmp,V+"Fig3_main.pdf",[(V+"Fig3D_heatmap.pdf",97,90)])
print(check(V+"Fig3_main.pdf",V+"Fig3_main_check.json"))
