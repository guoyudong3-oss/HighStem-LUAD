import sys; sys.path.insert(0,"figures/common")
from figlib import *
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/FigS4/"
fig=new_fig(190)
for t,x,y in [("a",1,1),("b",1,106)]: letter(fig,x,y,t)
tmp=V+"FigS4_v2_base.pdf"; fig.savefig(tmp,transparent=True)
embed_pdf(tmp,V+"FigS4_v2.pdf",[(V+"FigS4A_soft_power_v2.pdf",8,3),(V+"FigS4E_module_corr_v2.pdf",8,108)])
print(check(V+"FigS4_v2.pdf",V+"FigS4_v2_check.json"))
