# Shared helpers: 180 mm wide pages, text >= 8 pt, bold 10 pt panel letters, RGB without alpha
import os, re, subprocess, json
import matplotlib; matplotlib.use("Agg")
import matplotlib as mpl, matplotlib.pyplot as plt
import numpy as np
MM=1/25.4; W_MM=180.0
mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.labelsize":8,"axes.titlesize":8,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"legend.title_fontsize":8,
  "pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,"xtick.major.size":2.5,"ytick.major.size":2.5,"lines.linewidth":0.9,
  "figure.facecolor":"white","savefig.facecolor":"white","axes.unicode_minus":False})
def new_fig(h_mm):
    assert h_mm<=230.0,"height must not exceed 230 mm"
    f=plt.figure(figsize=(W_MM*MM,h_mm*MM)); f._h_mm=h_mm; return f
def ax_mm(fig,x,y,w,h,**kw):
    """x, y = top-left corner in mm from the page top-left; w, h in mm"""
    H=fig._h_mm; return fig.add_axes([x/W_MM,1-(y+h)/H,w/W_MM,h/H],**kw)
def letter(fig,x,y,t):
    """panel letter: bold 10 pt capital, top-left position in mm"""
    fig.text(x/W_MM,1-y/fig._h_mm,t.upper(),fontsize=10,fontweight="bold",ha="left",va="top",family="Arial")
def save(fig,out_base):
    fig.savefig(out_base+".pdf"); return out_base+".pdf"
def export_tiff(pdf,tiff,dpi):
    import fitz
    from PIL import Image
    Image.MAX_IMAGE_PIXELS=None
    pg=fitz.open(pdf)[0]; z=dpi/72; pix=pg.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False,colorspace=fitz.csRGB)
    Image.frombytes("RGB",(pix.width,pix.height),pix.samples).save(tiff,compression="tiff_lzw",dpi=(dpi,dpi)); return pix.width,pix.height
def render_png(pdf,png,dpi=200):
    import fitz
    pg=fitz.open(pdf)[0]; z=dpi/72; pg.get_pixmap(matrix=fitz.Matrix(z,z),alpha=False).save(png)
def check(pdf,report_json=None):
    """text checks: minimum font size, overlapping words, words touching the page edge, page size"""
    import fitz
    d=fitz.open(pdf); pg=d[0]; R=pg.rect; sizes=[]; small=[]; fonts=set()
    for b in pg.get_text("dict")["blocks"]:
        for l in b.get("lines",[]):
            for s in l["spans"]:
                if not s["text"].strip(): continue
                sizes.append(s["size"]); fonts.add(s["font"])
                if s["size"]<7.95: small.append((s["text"][:30],round(s["size"],2)))
    words=[w for w in pg.get_text("words") if w[4].strip()]
    over=[]
    for i in range(len(words)):
        a=fitz.Rect(words[i][:4])
        for j in range(i+1,len(words)):
            b=fitz.Rect(words[j][:4]); c=a&b
            if not c.is_empty and c.get_area()>0.3*min(a.get_area(),b.get_area()): over.append((words[i][4],words[j][4]))
    edge=[w[4] for w in words if w[0]<0.5 or w[1]<0.5 or w[2]>R.width-0.5 or w[3]>R.height-0.5]
    # cross-check with pdftotext -bbox (word-box height is about 1.1-1.2 x font size)
    hts=[]
    try:
        out=subprocess.run(["pdftotext","-bbox",pdf,"-"],capture_output=True,text=True).stdout
        for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>',out):
            wd=float(m.group(3))-float(m.group(1)); ht=float(m.group(4))-float(m.group(2))
            if wd>2.2*ht and len(m.group(5))>=4 and re.search(r"[A-Za-z0-9]",m.group(5)): hts.append(ht)   # horizontal alphanumeric words only
    except Exception as e: pass
    rep=dict(page_mm=[round(R.width/72*25.4,1),round(R.height/72*25.4,1)],min_font_pt=round(min(sizes),2) if sizes else None,max_font_pt=round(max(sizes),2) if sizes else None,
             n_text_spans=len(sizes),spans_below_8pt=small[:20],n_below_8=len(small),overlapping_word_pairs=over[:20],n_overlaps=len(over),words_touching_page_edge=edge[:20],
             pdftotext_min_word_height_pt=round(min(hts),2) if hts else None,fonts=sorted(fonts),images=len(pg.get_images()),drawings=len(pg.get_drawings()))
    if report_json: json.dump(rep,open(report_json,"w"),ensure_ascii=False,indent=1)
    return rep
def finalize(name,pdf,outdir,dpi=600,kind="mixed"):
    """copy the PDF, export TIFF (600 dpi mixed / 800 dpi line art) and a PNG preview, and write the check report"""
    import shutil
    os.makedirs(outdir,exist_ok=True); dst=os.path.join(outdir,name+".pdf")
    if os.path.abspath(pdf)!=os.path.abspath(dst): shutil.copy(pdf,dst)
    px=export_tiff(dst,os.path.join(outdir,name+".tif"),dpi); render_png(dst,os.path.join(outdir,name+"_render.png"),110)
    rep=check(dst,os.path.join(outdir,name+"_check.json")); rep.update(dpi=dpi,kind=kind,tif_pixels=list(px)); json.dump(rep,open(os.path.join(outdir,name+"_check.json"),"w"),ensure_ascii=False,indent=1); return rep

def embed_pdf(base_pdf,out_pdf,embeds):
    """embed existing PDFs (e.g. R heatmaps) 1:1 under the base page; embeds = [(src_pdf, x_mm, y_mm)]"""
    import fitz
    d=fitz.open(base_pdf); pg=d[0]
    for src,x,y in embeds:
        sd=fitz.open(src); sp=sd[0]; w,h=sp.rect.width,sp.rect.height; x0,y0=x/25.4*72,y/25.4*72
        pg.show_pdf_page(fitz.Rect(x0,y0,x0+w,y0+h),sd,0,overlay=False)
    d.save(out_pdf,deflate=True,garbage=3)
