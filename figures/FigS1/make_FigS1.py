import sys; sys.path.insert(0,"figures/common")
from figlib import *
V="results/reproduce_paper/runs/final_malignant_v3_20261009/W/FigS1/"
fig=new_fig(205); letter(fig,1,1,"a"); fig.savefig(V+"FigS1_v2_base.pdf",transparent=True)
embed_pdf(V+"FigS1_v2_base.pdf",V+"FigS1_v2.pdf",[(V+"out/FigS1.pdf",0,4.5)])
print(check(V+"FigS1_v2.pdf"))
