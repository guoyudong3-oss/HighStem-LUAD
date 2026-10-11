# Fig S7: median-split KM curves for the other five genes
import sys; sys.path.insert(0,"figures/Fig7")
from fig7lib import *
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/FigS8/"
fig=new_fig(136); rest=[g for g in GENES if g!="S100A16"]; res={}
for i,g in enumerate(rest):
    r,c=divmod(i,3); res[g]=km(fig,14+c*58,12+r*66,42,36,g,True,True)
for t,x,y in [("a",1,1),("b",59,1),("c",117,1),("d",1,67),("e",59,67)]: letter(fig,x,y,t)
out=V+"FigS8"; save(fig,out); print(res); print(check(out+".pdf",V+"FigS8_check.json"))
