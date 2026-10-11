# Collect the 8 main and 8 supplementary figures, export TIF and PNG previews, and check font size, clipping and overlap with poppler
import sys, os, re, shutil, subprocess, json
sys.path.insert(0,"figures/common")
import fitz
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
W="results/reproduce_paper/runs/final_malignant_v3_20261009/W/"; F=W+"FINAL/"
SRC={"Fig1":"Fig1/final/Fig1.pdf","Fig2":"Fig2/Fig2_main_v3.pdf","Fig3":"Fig3/final/Fig3.pdf","Fig4":"Fig4/final/Fig4.pdf","Fig5":"Fig5/final/Fig5.pdf","Fig6":"Fig6/Fig6_main_v3.pdf",
     "Fig7":"Fig7/final/Fig7.pdf","Fig8":"Fig8/Fig8_main_v2.pdf","FigS1":"FigS1/FigS1_v2.pdf","FigS2":"FigS2/FigS2_v2.pdf","FigS3":"FigS3/FigS3_v2.pdf","FigS4":"FigS4/FigS4_v2.pdf",
     "FigS5":"FigS5/FigS5_v2.pdf","FigS6":"FigS6/FigS6_v4.pdf","FigS7":"FigS8/final/FigS8.pdf","FigS8":"FigS9/final/FigS9.pdf"}
def poppler_check(pdf):
    xml=subprocess.run(["pdftohtml","-xml","-i","-stdout","-q","-zoom","3",pdf],capture_output=True,text=True).stdout
    specs={m.group(1):(int(m.group(2))/3.0,m.group(3)) for m in re.finditer(r'<fontspec id="(\d+)" size="(-?\d+)" family="([^"]+)"',xml)}
    used=set(re.findall(r'<text [^>]*font="(\d+)"',xml)); sizes=[specs[i][0] for i in used if i in specs]
    small=[(re.sub("<[^>]+>","",m.group(2))[:25],specs[m.group(1)][0]) for m in re.finditer(r'<text [^>]*font="(\d+)"[^>]*>(.*?)</text>',xml) if specs.get(m.group(1),(99,))[0]<7.9]
    bb=subprocess.run(["pdftotext","-bbox",pdf,"-"],capture_output=True,text=True).stdout
    pw,ph=[float(x) for x in re.search(r'<page width="([\d.]+)" height="([\d.]+)"',bb).groups()]
    words=[(float(a),float(b),float(c),float(d),t) for a,b,c,d,t in re.findall(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>',bb)]
    edge=[w[4] for w in words if w[0]<0.3 or w[1]<0.3 or w[2]>pw-0.3 or w[3]>ph-0.3]
    ov=[]
    for i in range(len(words)):
        a=words[i]
        for j in range(i+1,len(words)):
            b=words[j]; ix=max(0,min(a[2],b[2])-max(a[0],b[0])); iy=max(0,min(a[3],b[3])-max(a[1],b[1]))
            ar=lambda w:(w[2]-w[0])*(w[3]-w[1])
            if ix*iy>0.3*min(ar(a),ar(b)) and min(ar(a),ar(b))>0: ov.append((a[4],b[4]))
    hh=[w[3]-w[1] for w in words if (w[2]-w[0])>2.2*(w[3]-w[1]) and len(w[4])>=4]
    return dict(page_mm=(round(pw/72*25.4,1),round(ph/72*25.4,1)),min_font_pt=round(min(sizes),2) if sizes else None,below_8=small,edge_words=edge,overlaps=ov,min_word_h=round(min(hh),2) if hh else None,n_words=len(words),fonts=sorted({specs[i][1].split("+")[-1] for i in used if i in specs}))
res={}
for name,rel in SRC.items():
    dst=F+name+".pdf"; shutil.copy(W+rel,dst)
    d=fitz.open(dst); pg=d[0]; nimg=len(pg.get_images()); dpi=600 if nimg>0 else 800
    z=dpi/72; pix=pg.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False,colorspace=fitz.csRGB)
    Image.frombytes("RGB",(pix.width,pix.height),pix.samples).save(F+name+".tif",compression="tiff_lzw",dpi=(dpi,dpi))
    pg.get_pixmap(matrix=fitz.Matrix(150/72,150/72),alpha=False).save(F+"png_preview/"+name+".png")
    tif=Image.open(F+name+".tif"); txt=pg.get_text()
    r=poppler_check(dst); r.update(tif_dpi=[round(x) for x in tif.info.get("dpi",(0,0))],tif_mode=tif.mode,tif_px=tif.size,n_raster_images=nimg,kind="mixed" if nimg else "line art",old_numbers=[s for s in ("S7","S8","S9","Fig. S") if s in txt],source=rel)
    w_mm,h_mm=r["page_mm"]; r["pass"]=bool(r["min_font_pt"] is not None and r["min_font_pt"]>=7.9 and abs(w_mm-180)<=0.5 and h_mm<=230.5 and not r["edge_words"] and not r["overlaps"] and tif.mode=="RGB")
    res[name]=r; print(name,{k:r[k] for k in ("page_mm","min_font_pt","tif_dpi","tif_mode","edge_words","overlaps","pass")})
json.dump(res,open(F+"check_raw.json","w"),ensure_ascii=False,indent=1,default=str)
