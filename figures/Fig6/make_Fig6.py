# Fig6 A-E: single-gene ROC, repeated CV, development/test ROC and PR, confusion matrices, SHAP
import sys,glob,pickle; sys.path.insert(0,"figures/common"); sys.path.insert(0,"06_classifier")
from figlib import *
import pandas as pd, numpy as np
from sklearn.metrics import roc_curve, precision_recall_curve, auc
from fig6_model_labels import MODEL_LABELS
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/Fig6/"
cache=pickle.load(open(glob.glob(J+"J2/run_hgb/tables/fig6_model_cache_*.pkl")[0],"rb")); genes=[str(g) for g in cache["genes"]]; recs={r["name"]:r for r in cache["fit_records"]}
best=recs[cache["best_name"]]; ytr=np.asarray(cache["y_train"]).astype(int); yte=np.asarray(cache["y_test"]).astype(int); Xte=np.asarray(cache["X_test"],float); shap=np.asarray(cache["shap_vals"],float)
ptest=np.asarray(best["test_proba"],float); ptrain=np.asarray(best["train_proba"],float); correct=((ptest>=0.5).astype(int)==yte)
rc=np.load(J+"P/P3/panelC_roc_curves.npz"); tab10=plt.get_cmap("tab10")
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
fig=new_fig(228)
# A
a=ax_mm(fig,14,8,56,56)
for i,g in enumerate(genes):
    f,t=rc[g+"_fpr"],rc[g+"_tpr"]; a.plot(f,t,color=tab10(i),lw=1.0,label=f"{g}, AUC = {auc(f,t):.3f}")
a.plot([0,1],[0,1],"--",color="gray",lw=0.6); a.set_xlim(0,1); a.set_ylim(0,1.01); a.set_xlabel("1 - Specificity",fontsize=8,labelpad=1); a.set_ylabel("Sensitivity",fontsize=8,labelpad=4); clean(a); a.set_aspect("equal")
a.legend(loc="upper left",bbox_to_anchor=(-0.28,-0.2),ncol=2,columnspacing=0.8,frameon=False,fontsize=8,handlelength=1.4,labelspacing=0.25)
# B
pts=pd.read_csv(J+"J2/panel_D/Fig6_panel_D_benchmark_points.csv"); order=pts.groupby("model").auc.mean().sort_values().index.tolist()
b=ax_mm(fig,108,8,60,56); bp=b.boxplot([pts[pts.model==m].auc.values for m in order],vert=False,patch_artist=True,widths=0.6,showfliers=True,flierprops=dict(ms=1.5,mfc="none",mew=0.4),
 boxprops=dict(lw=0.5,facecolor="#9ECAE9"),medianprops=dict(color="black",lw=0.7),whiskerprops=dict(lw=0.5),capprops=dict(lw=0.5))
b.set_yticklabels([MODEL_LABELS[m]["display"] for m in order],fontsize=8); b.set_xlabel("AUC (5-fold CV x 10 repeats)",fontsize=8,labelpad=1); clean(b); b.set_xlim(0.95,1.0)
for i,m in enumerate(order,1): b.text(1.002,i,f"{pts[pts.model==m].auc.mean():.4f}",fontsize=8,va="center",ha="left",clip_on=False)
# C: development/test ROC and PR
col="#2AA7A7"
for r,(nm,y,p) in enumerate([("Development",ytr,ptrain),("Test",yte,ptest)]):
    fpr,tpr,_=roc_curve(y,p); pr,rc_,_=precision_recall_curve(y,p)
    a1=ax_mm(fig,14,98+r*37,36,22); a1.plot(fpr,tpr,color=col,lw=1.0); a1.plot([0,1],[0,1],"--",color="#BBBBBB",lw=0.6); a1.set_title(f"{nm}: AUC = {auc(fpr,tpr):.4f}",fontsize=8,pad=2,loc="left"); a1.set_xlabel("1 - Specificity",fontsize=8,labelpad=1); a1.set_ylabel("Sensitivity",fontsize=8,labelpad=1); clean(a1)
    a2=ax_mm(fig,68,98+r*37,30,22); a2.plot(rc_,pr,color=col,lw=1.0); a2.axhline(0.5,color="#BBBBBB",ls="--",lw=0.6); a2.set_title(f"PR AUC = {auc(rc_,pr):.4f}",fontsize=8,pad=2,loc="left"); a2.set_xlabel("Recall",fontsize=8,labelpad=1); a2.set_ylabel("Precision",fontsize=8,labelpad=1); clean(a2)
# D: confusion matrices
cf=pd.read_csv(J+"J2/table_confusion_HGB_forced.csv")
for k,(i,r) in enumerate(cf.iterrows()):
    ax=ax_mm(fig,124+k*31,104,22,22); M=np.array([[r.TN,r.FP],[r.FN,r.TP]]); tot=M.sum(1,keepdims=True)
    ax.imshow(M/tot,cmap="Blues",vmin=0,vmax=1.0,aspect="equal")
    for ii in range(2):
        for jj in range(2): ax.text(jj,ii,f"{M[ii,jj]:,}\n({M[ii,jj]/tot[ii,0]*100:.1f}%)",ha="center",va="center",fontsize=8,color="white" if M[ii,jj]/tot[ii,0]>0.55 else "black",linespacing=1.0)
    ax.set_xticks([0,1]); ax.set_xticklabels(["Low","High"],fontsize=8); ax.set_yticks([0,1]); ax.set_yticklabels(["Low","High"],fontsize=8,rotation=90,va="center")
    ax.set_xlabel("Predicted",fontsize=8,labelpad=1); (ax.set_ylabel("Observed",fontsize=8,labelpad=6) if k==0 else None)
    ax.set_title(("Development" if k==0 else "Test")+f"\naccuracy {r.accuracy:.4f}",fontsize=8,pad=2,linespacing=1.0)
    for sp in ax.spines.values(): sp.set_linewidth(0.4)
# E: SHAP
mean_abs=np.mean(np.abs(shap),axis=0); order=np.argsort(mean_abs); e=ax_mm(fig,36,182,86,38)
xmax=float(np.clip(np.quantile(np.abs(shap),0.992)*1.08,0.35,0.45)); rng=np.random.default_rng(42); cmap=plt.get_cmap("rainbow"); norm=mpl.colors.PowerNorm(gamma=1.35,vmin=0,vmax=1)
for pos,idx in enumerate(order):
    x=shap[:,idx]; fv=Xte[:,idx]; c=(fv-fv.min())/(fv.max()-fv.min()); y=pos+rng.normal(0,0.085,len(x))
    e.scatter(x[correct],y[correct],c=c[correct],cmap=cmap,norm=norm,s=2.5,marker="o",alpha=0.7,linewidths=0,rasterized=True)
    e.scatter(x[~correct],y[~correct],c=c[~correct],cmap=cmap,norm=norm,s=5,marker="x",alpha=0.6,linewidths=0.3,rasterized=True)
    e.text(-xmax+0.01*2*xmax,pos+0.38,f"{mean_abs[idx]:.3f}",fontsize=8,color="#333333",va="center")
e.axvline(0,color="#8A8A8A",lw=0.6); e.set_xlim(-xmax,xmax); e.set_yticks(range(len(genes))); e.set_yticklabels([genes[i] for i in order],fontsize=8); e.set_xlabel("SHAP value (impact on HighStem probability)",fontsize=8,labelpad=1); clean(e)
cax=ax_mm(fig,128,188,2.4,22); cb=fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap,norm=norm),cax=cax,ticks=[0,1]); cb.ax.set_yticklabels(["Low","High"]); cb.ax.tick_params(labelsize=8,length=2,width=0.5); cb.outline.set_linewidth(0.4)
fig.text(127/W_MM,1-185/fig._h_mm,"Expression",fontsize=8,va="center")
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],color="black",marker="x",ls="None",ms=3.5,label="Misclassified"),Line2D([],[],color="black",marker="o",ls="None",ms=3.5,label="Correct")],loc="upper left",bbox_to_anchor=(138/W_MM,1-188/fig._h_mm),frameon=False,fontsize=8,borderaxespad=0)
fig.text(36/W_MM,1-179/fig._h_mm,"Numbers: mean |SHAP|",fontsize=8,va="center")
for t,x,y in [("a",1,1),("b",92,1),("c",1,92),("d",108,92),("e",1,174)]: letter(fig,x,y,t)
out=V+"Fig6_main_v3"; save(fig,out)
