# Fig5 A-E: activity genes, spatial projection (embedded R output), TCGA tumour vs normal, KM (optimal cutoff), UpSet
import sys; sys.path.insert(0,"figures/common")
from figlib import *
import pandas as pd, numpy as np
from itertools import combinations
from scipy.stats import gaussian_kde
from lifelines import KaplanMeierFitter
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/Fig5/"; T=J+"P/P2/run_fig5/tables/"
fig=new_fig(206)
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
# A: genes ranked by Pearson r (top 20 labelled)
d=pd.read_csv(J+"O/O3/Fig5A_all138_genes_r_p.csv").sort_values("rank_by_r"); n=len(d)
a=ax_mm(fig,14,10,62,96); y=-np.arange(n); lp=-np.log10(d.p.values)
cm=mpl.colors.LinearSegmentedColormap.from_list("pb",["#9ECAE9","#2C7FB8","#132B43"])
a.hlines(y,0,d.r,color="#E3E3E3",lw=0.3,zorder=1)
sc=a.scatter(d.r,y,c=lp,cmap=cm,s=7,linewidths=0,vmin=lp.min(),vmax=lp.max(),zorder=3)
a.set_ylim(-n-1,3); a.set_xlim(0.14,0.5); a.set_yticks([]); a.set_xlabel("Pearson r (HighStem cells)",fontsize=8,labelpad=1); a.set_xticks([0.15,0.2,0.25,0.3]); clean(a); a.spines["left"].set_visible(False)
top=d.head(20); ly=np.linspace(1,-98,20)
for (i,r),yy in zip(top.iterrows(),ly):
    a.text(0.345,yy,r.gene,fontsize=8,va="center",ha="left"); a.plot([r.r+0.002,0.341],[-int(r.rank_by_r)+1,yy],color="#999999",lw=0.3,zorder=2)
cax=ax_mm(fig,14,121,40,2.6); cb=fig.colorbar(sc,cax=cax,orientation="horizontal"); cb.set_ticks([20,40,60,80]); cb.set_ticklabels(["1e-20","1e-40","1e-60","1e-80"]); cb.ax.tick_params(labelsize=8,length=2,width=0.5,pad=1); cb.outline.set_linewidth(0.4)
fig.text(58/W_MM,1-122.3/fig._h_mm,"p value",fontsize=8,va="bottom")
fig.text(14/W_MM,1-5/fig._h_mm,"138 activity genes (top 20 labelled)",fontsize=8,fontweight="bold",va="center")
# C: TCGA tumour vs normal
sc_all=pd.read_csv(T+"tcga_highstem_activity_scores_all_fromA_primary517.csv"); cc=ax_mm(fig,118,66,38,40)
cols={"Normal":"#2A9D9D","Tumor":"#D6204A"}; rng=np.random.default_rng(2)
for i,g in enumerate(["Normal","Tumor"]):
    v=sc_all[sc_all.sample_type==g].PlasticStem_Score.values; kd=gaussian_kde(v); ys=np.linspace(v.min(),v.max(),200); dd=kd(ys); dd=dd/dd.max()*0.42
    cc.fill_betweenx(ys,i-dd,i+dd,color=cols[g],alpha=0.45,lw=0)
    cc.scatter(i+rng.uniform(-0.12,0.12,len(v)),v,s=0.8,color=cols[g],linewidths=0,rasterized=True)
    cc.boxplot([v],positions=[i],widths=0.12,showfliers=False,patch_artist=True,boxprops=dict(facecolor="white",lw=0.5),medianprops=dict(color="black",lw=0.7),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5))
yt=sc_all.PlasticStem_Score.max()+0.1; cc.plot([0,0,1,1],[yt-0.05,yt,yt,yt-0.05],color="black",lw=0.6); cc.text(0.5,yt+0.02,"****",fontsize=8,ha="center",va="bottom")
cc.set_xticks([0,1]); cc.set_xticklabels([f"Normal\n(n = 59)","Tumor\n(n = 517)"]); cc.set_ylabel("HighStem activity score",fontsize=8,labelpad=1); cc.set_xlim(-0.6,1.6); cc.set_ylim(sc_all.PlasticStem_Score.min()-0.1,yt+0.45); clean(cc)
# D: KM curves and risk table
sv=pd.read_csv(T+"tcga_highstem_activity_scores_tumor01_fromA.csv").dropna(subset=["OS","OS.time"]); cut=0.663686886636902; sv["t"]=sv["OS.time"]/365.25; sv["hi"]=sv.PlasticStem_Score>cut+1e-9
assert int(sv.hi.sum())==78 and int((~sv.hi).sum())==430,(sv.hi.sum(),)
k=ax_mm(fig,14,146,76,34); 
for g,c,lab in [(True,"#C0392B","High"),(False,"#1B7F5C","Low")]:
    m=sv[sv.hi==g]; kmf=KaplanMeierFitter().fit(m.t,m["OS"]); kmf.plot_survival_function(ax=k,ci_show=False,color=c,lw=1.0,show_censors=True,censor_styles=dict(ms=3,marker="|",mew=0.6),label=lab)
k.set_xlim(0,10); k.set_ylim(0,1.02); k.set_xticks([0,2.5,5,7.5,10]); k.set_xticklabels(["0","2.5","5","7.5","10"]); k.set_yticks([0,0.25,0.5,0.75,1.0]); k.set_yticklabels(["0","0.25","0.5","0.75","1.0"])
k.set_xlabel("Time (years)",fontsize=8,labelpad=1); k.set_ylabel("Survival probability",fontsize=8,labelpad=1); k.legend(frameon=False,fontsize=8,loc="upper right",handlelength=1.2,title="Score",title_fontsize=8); clean(k)
k.text(0.97,0.58,"Log-rank p = 0.0068",transform=k.transAxes,fontsize=8,va="center",ha="right")
tb=ax_mm(fig,14,194,76,8.5); times=[0,2.5,5,7.5,10]
for yi,(g,lab,c) in enumerate([(True,"High","#C0392B"),(False,"Low","#1B7F5C")]):
    m=sv[sv.hi==g]
    for tt in times: tb.text(tt,1-yi,str(int((m.t>=tt).sum())),fontsize=8,ha="center",va="center")
    tb.text(-0.4,1-yi,lab,fontsize=8,ha="right",va="center",color=c)
tb.set_xlim(0,10); tb.set_ylim(-0.6,1.6); tb.axis("off"); fig.text(14/W_MM,1-190.5/fig._h_mm,"Number at risk",fontsize=8,va="bottom") if False else None
tb2=ax_mm(fig,14,189.5,76,3); tb2.axis("off"); tb2.text(0,0.5,"Number at risk",fontsize=8,va="center")
# E: UpSet
sets={"Boruta":"ml_selected_genes_Boruta.csv","LASSO":"ml_selected_genes_LASSO.csv","ABESS":"ml_selected_genes_ABESS_Proxy.csv","GBM":"ml_selected_genes_GBM.csv","DT":"ml_selected_genes_DT.csv","RF":"ml_selected_genes_RF.csv"}
S={k:set(pd.read_csv(T+v).gene) for k,v in sets.items()}; names=list(S); allg=set().union(*S.values())
from collections import Counter
mem={g:tuple(int(g in S[k]) for k in names) for g in allg}; cnt=Counter(mem.values()); items=cnt.most_common(15)
items=sorted(items,key=lambda x:-x[1])
ub=ax_mm(fig,127,146,50,24); ub.bar(range(len(items)),[c for _,c in items],color="#222222",width=0.7)
for i,(_,c) in enumerate(items): ub.text(i,c+1,str(c),fontsize=8,ha="center",va="bottom")
ub.set_xlim(-0.6,len(items)-0.4); ub.set_ylim(0,max(c for _,c in items)*1.2); ub.set_xticks([]); ub.set_ylabel("Intersection size",fontsize=8,labelpad=1); clean(ub)
um=ax_mm(fig,127,172,50,24); um.set_xlim(-0.6,len(items)-0.4); um.set_ylim(-0.6,len(names)-0.4); um.axis("off")
for yi in range(len(names)): um.add_patch(plt.Rectangle((-0.6,yi-0.5),len(items),1,color="#F4F4F4" if yi%2==0 else "white",lw=0,zorder=0))
for i,(m,_) in enumerate(items):
    ys=[len(names)-1-j for j,b in enumerate(m) if b]
    for j,b in enumerate(m): um.scatter(i,len(names)-1-j,s=9,color="black" if b else "#CCCCCC",zorder=3,linewidths=0)
    if len(ys)>1: um.plot([i,i],[min(ys),max(ys)],color="black",lw=0.8,zorder=2)
sb=ax_mm(fig,94,172,16,24); sb.barh([len(names)-1-i for i in range(len(names))],[len(S[k]) for k in names],color="#555555",height=0.6); sb.invert_xaxis(); sb.set_ylim(-0.6,len(names)-0.4); sb.set_yticks([]); sb.set_xticks([0,60,120]); sb.set_xlabel("Set size",fontsize=8,labelpad=1); clean(sb); sb.spines["left"].set_visible(False)
for j,k in enumerate(names): fig.text(126/W_MM,1-(172+(j+0.5)*24/len(names))/fig._h_mm,k,fontsize=8,ha="right",va="center")
six=sorted(g for g,m in mem.items() if all(m)); i6=[i for i,(m,_) in enumerate(items) if all(m)][0]
fig.text(158/W_MM,1-154/fig._h_mm,"\n".join(six),fontsize=8,va="center",ha="left",linespacing=1.05)
for t,x,y in [("a",1,1),("b",92,1),("c",110,60),("d",1,140),("e",92,140)]: letter(fig,x,y,t)
tmp=V+"Fig5_base.pdf"; fig.savefig(tmp,transparent=True)
embed_pdf(tmp,V+"Fig5_main.pdf",[(V+"panels/Fig5B_spatial.pdf",92,6)])
print(check(V+"Fig5_main.pdf",V+"Fig5_main_check.json"))
