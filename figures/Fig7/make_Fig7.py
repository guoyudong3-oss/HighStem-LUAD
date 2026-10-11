# Fig7: A tumour vs normal, B paired samples, C ROC, D S100A16 KM (TCGA-LUAD: 517 tumours, 59 normals, 58 pairs, 508 with survival)
import sys; sys.path.insert(0,"figures/Fig7")
from fig7lib import *
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/Fig7/"
fig=new_fig(186); PW=22; X0=[14+i*28.8 for i in range(6)]
for i,g in enumerate(GENES):
    a=ax_mm(fig,X0[i],10,PW,36); nrm=expr_use.loc[expr_use.type=="Normal",g].astype(float).values; tu=expr_use.loc[expr_use.type=="Tumor",g].astype(float).values
    parts=a.violinplot([nrm,tu],positions=[0,1],widths=0.85,showextrema=False)
    for j,bd in enumerate(parts["bodies"]): bd.set_alpha(0.38); bd.set_facecolor("#149C9A" if j==0 else "#D87D4A"); bd.set_edgecolor("none")
    bp=a.boxplot([nrm,tu],positions=[0,1],widths=0.16,patch_artist=True,showfliers=False,boxprops=dict(lw=0.5),medianprops=dict(color="black",lw=0.7),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5)); bp["boxes"][0].set_facecolor("#BFDCDC"); bp["boxes"][1].set_facecolor("#F2B38E")
    rng=np.random.default_rng(42); a.scatter(rng.normal(0,0.05,len(nrm)),nrm,s=1.5,color="#108A88",linewidths=0,rasterized=True); a.scatter(rng.normal(1,0.05,len(tu)),tu,s=1.5,color="#C81D25",linewidths=0,rasterized=True)
    p=mannwhitneyu(nrm,tu,alternative="two-sided").pvalue; ymin=min(nrm.min(),tu.min()); ymax=max(nrm.max(),tu.max()); sp=ymax-ymin
    a.set_ylim(ymin-0.05*sp,ymax+0.22*sp); a.text(0.5,ymax+0.06*sp,fp(p),fontsize=8,ha="center",va="bottom"); a.set_xticks([0,1]); a.set_xticklabels(["Normal","Tumor"],rotation=40,ha="right",rotation_mode="anchor"); a.set_title(g,fontsize=8,fontweight="bold",pad=2)
    if i==0: a.set_ylabel("log2(TPM + 1)",fontsize=8,labelpad=1)
    clean(a)
    b=ax_mm(fig,X0[i],66,PW,36); pv=expr_use.pivot_table(index="patient",columns="type",values=g,aggfunc="mean").dropna()
    for _,r in pv.iterrows(): b.plot([0,1],[r.Normal,r.Tumor],color="black",alpha=0.6,lw=0.5)
    bp=b.boxplot([pv.Normal.values,pv.Tumor.values],positions=[0,1],widths=0.4,patch_artist=True,showfliers=False,boxprops=dict(lw=0.5),medianprops=dict(color="black",lw=0.7),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5)); bp["boxes"][0].set_facecolor("#149C9A"); bp["boxes"][1].set_facecolor("#C81D25")
    pp=wilcoxon(pv.Normal,pv.Tumor).pvalue; ymin=pv.values.min(); ymax=pv.values.max(); sp=ymax-ymin; b.set_ylim(ymin-0.05*sp,ymax+0.22*sp); b.text(0.5,ymax+0.06*sp,fp(pp),fontsize=8,ha="center",va="bottom")
    b.set_xticks([0,1]); b.set_xticklabels(["Normal","Tumor"],rotation=40,ha="right",rotation_mode="anchor"); b.set_title(g,fontsize=8,fontweight="bold",pad=2)
    if i==0: b.set_ylabel("log2(TPM + 1)",fontsize=8,labelpad=1)
    clean(b)
# C ROC
c=ax_mm(fig,14,122,50,38); y=(expr_use.type=="Tumor").astype(int).values; pal=["#273C75","#C0392B","#D4A62A","#2E8B57","#6AB7BE","#8E44AD"]
for i,g in enumerate(GENES):
    s=expr_use[g].astype(float).values; f,t,_=roc_curve(y,s); c.plot(f,t,lw=1.0,color=pal[i],label=f"{g}, AUC = {roc_auc_score(y,s):.3f}")
c.plot([0,1],[0,1],color="#BBBBBB",lw=0.6); c.set_xlim(0,1); c.set_ylim(0,1.02); c.set_xlabel("1 - Specificity",fontsize=8,labelpad=1); c.set_ylabel("Sensitivity",fontsize=8,labelpad=1); clean(c); c.legend(loc="upper left",bbox_to_anchor=(-0.1,-0.22),ncol=2,frameon=False,fontsize=8,handlelength=1.2,labelspacing=0.25,columnspacing=0.8)
# D KM S100A16
km(fig,108,122,56,38,"S100A16",True,True)
for t,x,y_ in [("a",1,1),("b",1,57),("c",1,113),("d",90,113)]: letter(fig,x,y_,t)
out=V+"Fig7_main"; save(fig,out); print(check(out+".pdf",V+"Fig7_main_check.json"))
