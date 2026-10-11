# Fig4 A-D: dendrogram, kME, module UMAP, module dot plot (R panels embedded at final size)
import sys; sys.path.insert(0,"figures/common")
from figlib import *
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/Fig4/"; P=V+"panels/"
fig=new_fig(222)
for t,x,y in [("a",1,1),("b",1,64),("c",1,129),("d",124,129)]: letter(fig,x,y,t)
tmp=V+"Fig4_base.pdf"; fig.savefig(tmp,transparent=True)
embed_pdf(tmp,V+"Fig4_main.pdf",[(P+"Fig4A_dendrogram.pdf",0,3),(P+"Fig4B_kME.pdf",0,66),(P+"Fig4C_module_umap.pdf",0,132),(P+"Fig4D_dotplot.pdf",122,132)])
print(check(V+"Fig4_main.pdf",V+"Fig4_main_check.json"))
