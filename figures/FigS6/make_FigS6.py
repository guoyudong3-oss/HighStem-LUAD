# Fig S6: density maps, model PR/ROC, decision curve, dependence plots
import sys,glob,pickle; sys.path.insert(0,"figures/common"); sys.path.insert(0,"06_classifier"); sys.path.insert(0,"06_classifier")
from figlib import *
import pandas as pd, numpy as np
from scipy.stats import gaussian_kde
from matplotlib.path import Path
from sklearn.metrics import roc_curve, precision_recall_curve
from fig6_model_labels import MODEL_LABELS
import make_fig6_panel_ef_only as EF
J="results/reproduce_paper/runs/final_malignant_v3_20261009/"; V=J+"W/FigS6/"
d=pickle.load(open(V+"data/ab_data.pkl","rb")); genes=d["genes"]
cache=pickle.load(open(glob.glob(J+"J2/run_hgb/tables/fig6_model_cache_*.pkl")[0],"rb")); recs={r["name"]:r for r in cache["fit_records"]}; best=recs[cache["best_name"]]
yte=np.asarray(cache["y_test"]).astype(int); Xte=np.asarray(cache["X_test"],float); shap=np.asarray(cache["shap_vals"],float); ptest=np.asarray(best["test_proba"],float); correct=((ptest>=0.5).astype(int)==yte)
def clean(a): a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
fig=new_fig(230)
def dens(ax,xy,w,xlim,ylim,outline,title):
    x0,x1=xlim; y0,y1=ylim; xx,yy=np.mgrid[x0:x1:160j,y0:y1:160j]; grid=np.vstack([xx.ravel(),yy.ravel()]); wp=np.maximum(w,0)
    if float(wp.sum())<=0: wp=np.ones_like(wp)
    zz=gaussian_kde(xy.T,weights=wp).__call__(grid).reshape(xx.shape)   # all cells, as in weighted_density of the original script
    polys=outline if isinstance(outline,list) else [outline]
    mask=np.zeros(xx.shape,bool)
    for p in polys:
        if p is not None and p.shape[0]>2: mask|=Path(p).contains_points(np.c_[xx.ravel(),yy.ravel()]).reshape(xx.shape)
    zz=np.ma.array(zz,mask=~mask) if mask.any() else zz
    m=ax.contourf(xx,yy,zz,levels=18,cmap="magma"); ax.set_facecolor("black")
    for p in polys:
        if p is not None and p.shape[0]>2: ax.plot(p[:,0],p[:,1],color="white",ls="--",lw=0.6)
    ax.set_xlim(x0,x1); ax.set_ylim(y0,y1); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title,fontsize=8,fontweight="bold",pad=1.5)
    return m
for row,(coords,expr,xl,yl,ol,y) in enumerate([(d["coords_a"],d["expr"],d["xlim_a"],d["ylim_a"],d["outline_a"],6),(d["coords_b"],d["expr_malignant"],d["xlim_b"],d["ylim_b"],d["outline_b"],50)]):
    for i,g in enumerate(genes):
        ax=ax_mm(fig,6+i*28.5,y+5,18,34); m=dens(ax,coords,expr[:,i],xl,yl,ol,g)
        cb=fig.colorbar(m,ax=ax,fraction=0.07,pad=0.03,aspect=18); cb.ax.tick_params(labelsize=8,length=1.5,width=0.4); cb.outline.set_linewidth(0.3); vmx=float(np.nanmax(m.levels)); ex=10**np.floor(np.log10(vmx)); top=np.floor(vmx/ex)*ex; cb.set_ticks([0,top]); cb.set_ticklabels(['0',f'{top:.1g}'])
# model PR/ROC
pal=EF.ef_model_palette(); a1=ax_mm(fig,14,108,30,28); a2=ax_mm(fig,56,108,30,28)
for k in EF.ef_model_order():
    p=np.asarray(recs[k]["test_proba"],float); pr,rc,_=precision_recall_curve(yte,p); fpr,tpr,_=roc_curve(yte,p)
    a1.plot(rc,pr,color=pal[k],lw=0.9,label=MODEL_LABELS[k]["display"]); a2.plot(fpr,tpr,color=pal[k],lw=0.9)
a1.set_xlabel("Recall",fontsize=8,labelpad=1); a1.set_ylabel("Precision",fontsize=8,labelpad=1); a2.set_xlabel("1 - Specificity",fontsize=8,labelpad=1); a2.set_ylabel("Sensitivity",fontsize=8,labelpad=1)
a1.set_ylim(0.9,1.005); a1.set_xlim(0,1.0); a2.set_xlim(0,1); a2.set_ylim(0,1.01); a2.plot([0,1],[0,1],"--",color="#BBBBBB",lw=0.5); clean(a1); clean(a2)
fig.legend(*a1.get_legend_handles_labels(),loc="center left",bbox_to_anchor=(90/W_MM,1-122/fig._h_mm),frameon=False,fontsize=8,handlelength=1.3,labelspacing=0.3,title="Model",title_fontsize=8)
# H DCA
dca=pd.read_csv(J+"J2/table_dca_HGB_forced.csv"); h=ax_mm(fig,150,108,27,28)

h.plot(dca.threshold,dca.net_benefit_model,color="#D62728",lw=1.0,label="Model"); h.plot(dca.threshold,dca.net_benefit_treat_all,color="#7F7F7F",lw=0.8,label="Treat all"); h.plot(dca.threshold,dca.net_benefit_treat_none,color="black",lw=0.8,ls="--",label="Treat none")
h.set_xlabel("Threshold probability",fontsize=8,labelpad=1); h.set_ylabel("Net benefit",fontsize=8,labelpad=1); h.set_ylim(-0.05,0.55); h.legend(frameon=False,fontsize=8,loc="center right",handlelength=1.2); clean(h)
# normalised expression vs predicted probability / SHAP
def dep(x0,y0,vals,ylab,ttl):
    for i,g in enumerate(genes):
        r,c=divmod(i,3); ax=ax_mm(fig,x0+c*26,y0+r*32,19,21); x=Xte[:,i]; x=(x-x.min())/(x.max()-x.min()); v=vals(i)
        ax.scatter(x[correct],v[correct],s=1.2,color="#555555",alpha=0.5,linewidths=0,rasterized=True); ax.scatter(x[~correct],v[~correct],s=4,marker="^",color="black",alpha=0.7,linewidths=0,rasterized=True)
        ax.set_title(g,fontsize=8,pad=1.5); ax.set_xticks([0,0.5,1]); ax.set_xticklabels(["0","0.5","1"])
        if c==0: ax.set_ylabel(ylab,fontsize=8,labelpad=1)
        if r==1: ax.set_xlabel("Expression (0-1)",fontsize=8,labelpad=1)
        clean(ax)
dep(14,160,lambda i:ptest,"P(HighStem)","J"); dep(104,160,lambda i:shap[:,i],"SHAP value","K")
fig.legend(handles=[mpl.lines.Line2D([],[],color="#555555",marker="o",ls="None",ms=3,label="Correct"),mpl.lines.Line2D([],[],color="black",marker="^",ls="None",ms=3.5,label="Misclassified")],loc="lower center",ncol=2,bbox_to_anchor=(0.5,0.0),frameon=False,fontsize=8)
for t,x,y in [("a",1,1),("b",1,45),("c",1,100),("d",142,100),("e",1,152),("f",92,152)]: letter(fig,x,y,t)
out=V+"FigS6_v4"; save(fig,out)
