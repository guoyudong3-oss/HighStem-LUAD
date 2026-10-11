import pandas as pd
from collections import Counter
R="results/reproduce_paper/runs/"; N=R+"final_malignant_v2_20261009/I1/run/tables/"; L=R+"paper_main_epiref_latest13k_originalflow_20260214_125353/tables/"; R3=R+"rerun3_latest13k_20261009/tables/"
act=pd.read_csv(N+"high_stem_activity_genes_final.csv"); print("Pearson passed",len(act),"r range",round(act.r.min(),3),round(act.r.max(),3),"max p",act.p.max())
print("top 10:",act.head(10)[["gene","r"]].round(3).values.tolist())
algs=["Boruta","RF","DT","LASSO","GBM","ABESS_Proxy"]
sel={a:pd.read_csv(N+f"ml_selected_genes_{a}.csv").gene.tolist() for a in algs}
lat={a:set(pd.read_csv(L+f"ml_selected_genes_{a}.csv").gene) for a in algs}
latact=pd.read_csv(L+"high_stem_activity_genes_final.csv")
print("same activity gene set as the originalflow run:",set(act.gene)==set(latact.gene),len(latact))
for a in algs: print(f"{a:12s} n={len(sel[a]):3d} same list as originalflow:",set(sel[a])==lat[a],"(originalflow n=%d)"%len(lat[a]))
c=Counter(g for v in sel.values() for g in v); vote=pd.DataFrame(sorted(c.items(),key=lambda x:(-x[1],x[0])),columns=["gene","votes"])
for a in algs: vote[a]=vote.gene.isin(sel[a]).astype(int)
vote.to_csv(N.replace("run/tables/","")+"vote_table_with_membership.csv",index=False)
print("vote distribution",vote.votes.value_counts().sort_index(ascending=False).to_dict())
print(vote[vote.votes>=4].to_string(index=False))
six=vote.gene[vote.votes==6].tolist(); r3=pd.read_csv(R3+"final_consensus_genes.csv").gene.tolist(); old5=["S100A16","MARCKSL1","MUC5B","TIMP1","CST6"]
print("6 votes",six); print("overlap with rerun3 five genes",sorted(set(six)&set(r3)),"| overlap with earlier five genes",sorted(set(six)&set(old5)))
print("votes of rerun3 five genes",{g:int(c.get(g,0)) for g in r3},"| votes of earlier five genes",{g:int(c.get(g,0)) for g in old5})
print("originalflow final_consensus:",pd.read_csv(L+"final_consensus_genes.csv").gene.tolist(),"| this run:",pd.read_csv(N+"final_consensus_genes.csv").gene.tolist())
for a in algs: print(a,sorted(sel[a]))
