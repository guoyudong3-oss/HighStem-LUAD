# Export Fig6 tables from the model cache (threshold 0.5; DCA thresholds linspace(0.01, 0.99, 181)); labels use sklearn class names
import pickle,glob,sys,numpy as np,pandas as pd
sys.path.insert(0,"06_classifier"); sys.path.insert(0,"06_classifier")
from make_fig6_reference_like import net_benefit
from fig6_model_labels import MODEL_LABELS
from sklearn.metrics import confusion_matrix,roc_auc_score,average_precision_score,precision_recall_curve,auc
J="results/reproduce_paper/runs/final_malignant_v3_20261009/J2/"
for tag,run in [("HGB_forced","run_hgb"),("RF_best_by_CV","run")]:
    c=pickle.load(open(glob.glob(J+run+"/tables/fig6_model_cache_*.pkl")[0],"rb")); genes=[str(g) for g in c["genes"]]
    ytr=np.asarray(c["y_train"]).astype(int); yte=np.asarray(c["y_test"]).astype(int)
    recs={r["name"]:r for r in c["fit_records"]}
    rows=[]
    for k,r in recs.items():
        row=dict(model=k,model_display=MODEL_LABELS[k]["display"],sklearn_class=MODEL_LABELS[k]["sklearn_class"])
        for sp,y,p in [("dev",ytr,np.asarray(r["train_proba"])),("test",yte,np.asarray(r["test_proba"]))]:
            fpr_p,rec_p,_=precision_recall_curve(y,p); row[f"{sp}_ROC_AUC"]=roc_auc_score(y,p); row[f"{sp}_PR_AUC"]=auc(rec_p,fpr_p)
        row["cv_mean_auc_inner5fold"]=float(np.mean(r["cv_auc"])); rows.append(row)
    pd.DataFrame(rows).to_csv(J+f"table_all_models_dev_test_AUC_{tag}.csv",index=False)
    best=c["best_name"]; r=recs[best]; ptr=np.asarray(r["train_proba"]); pte=np.asarray(r["test_proba"])
    cm=[]
    for sp,y,p in [("development(train)",ytr,ptr),("test",yte,pte)]:
        tn,fp,fn,tp=confusion_matrix(y,(p>=0.5).astype(int),labels=[0,1]).ravel(); cm.append(dict(split=sp,n=len(y),TP=tp,FP=fp,FN=fn,TN=tn,accuracy=(tp+tn)/len(y),model=MODEL_LABELS[best]["display"]))
    pd.DataFrame(cm).to_csv(J+f"table_confusion_{tag}.csv",index=False)
    th=np.linspace(0.01,0.99,181); prev=yte.mean()
    pd.DataFrame(dict(threshold=th,net_benefit_model=net_benefit(yte,pte,th),net_benefit_treat_all=prev-(1-prev)*(th/(1-th)),net_benefit_treat_none=0.0)).to_csv(J+f"table_dca_{tag}.csv",index=False)
    sh=np.asarray(c["shap_vals"],float); s=pd.DataFrame(dict(gene=genes,mean_abs_shap=np.abs(sh).mean(0),mean_shap=sh.mean(0),model=MODEL_LABELS[best]["display"])).sort_values("mean_abs_shap",ascending=False); s.to_csv(J+f"table_shap_{tag}.csv",index=False)
    print("==",tag,"model",MODEL_LABELS[best]["display"],"genes",genes); print(pd.DataFrame(rows)[["model_display","dev_ROC_AUC","dev_PR_AUC","test_ROC_AUC","test_PR_AUC","cv_mean_auc_inner5fold"]].round(4).to_string(index=False)); print(pd.DataFrame(cm).round(4).to_string(index=False)); print(s.round(4).to_string(index=False)); print("DCA test max NB",np.nanmax(net_benefit(yte,pte,th)).round(4),"prevalence",round(prev,4))
