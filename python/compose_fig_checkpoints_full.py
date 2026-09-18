#!/usr/bin/env python3
"""Combined PD-L1 / PD-1 checkpoint figure: PD-L1 immunohistochemistry (images + quantification)
together with the flow-cytometry checkpoint data (PD-L1 and PD-1 gates + time courses).
  A  PD-L1 IHC images            B  PD-L1 IHC quantification (48 h)
  C  PD-L1 flow gate | PD-L1+ (% of CD45- spheroid) over time
  D  PD-1  flow gate | PD-1+  (% of CD3+ T cells) over time
Output: output/figuras/Fig_checkpoints_full.png"""
import os, numpy as np
from PIL import Image, ImageDraw, ImageFont

LAY = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/LAYOUTS VIABILIDAD MANUSCRITO"
FIG = "output/figuras"
FDIR = "/usr/share/fonts/truetype/dejavu"
def font(sz, bold=True):
    return ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"), sz)

LAYS = {
    "PDL1": ("PD-L1", "PD-L1+ cells", "CD45- spheroid", "FSC-A", "PD-L1"),
    "PD1":  ("PD-1+ ", "PD-1+ T cells", "CD3+ (CD45+)", "CD3", "PD-1"),
}

def load_plot(name):
    im = Image.open(os.path.join(LAY, name + ".png")).convert("RGBA")
    im = Image.alpha_composite(Image.new("RGBA", im.size, (255, 255, 255, 255)), im).convert("RGB")
    g = np.array(im.convert("L"))
    r = np.where(g.max(axis=1) > 150)[0]; c = np.where(g.max(axis=0) > 150)[0]
    im = im.crop((c.min(), r.min(), c.max() + 1, r.max() + 1)); g = g[r.min():r.max()+1, c.min():c.max()+1]
    H, Wd = g.shape
    vc = [i for i in range(Wd) if (g[:, i] < 80).sum() > 0.6 * H]
    hr = [i for i in range(H) if (g[i, :] < 80).sum() > 0.6 * Wd]
    return im.crop((min(vc), min(hr), max(vc) + 1, max(hr) + 1))

plots = {k: load_plot(v[0]) for k, v in LAYS.items()}
# relabel gate stats -> English + period decimals (as in compose_fig_checkpoints)
_pl = plots["PDL1"]; _dp = ImageDraw.Draw(_pl); _pf = font(28, bold=False)
_dp.rectangle([54, 168, 172, 242], fill="white")
_dp.text((168, 176), "PD-L1+", fill="black", font=_pf, anchor="ra")
_dp.text((168, 208), "9.58", fill="black", font=_pf, anchor="ra")
_p1 = plots["PD1"]; _d1 = ImageDraw.Draw(_p1); _1f = font(28, bold=False)
_d1.rectangle([278, 158, 434, 230], fill="white")
_d1.text((282, 166), "PD-1+", fill="black", font=_1f, anchor="la")
_d1.text((282, 198), "89.7", fill="black", font=_1f, anchor="la")

pw = max(p.width for p in plots.values()); ph = max(p.height for p in plots.values())
ylab_w, header_h, xlab_h = 52, 74, 44
cell_w, cell_h = ylab_w + pw, header_h + ph + xlab_h
fhead, fpar, fax = font(30), font(23, False), font(28)

def layout_cell(key):
    _, title, parent, xl, yl = LAYS[key]
    cell = Image.new("RGB", (cell_w, cell_h), "white"); dd = ImageDraw.Draw(cell)
    dd.text((2, 0), title, fill="black", font=fhead, anchor="la")
    dd.text((2, 38), "from " + parent, fill=(90, 90, 90), font=fpar, anchor="la")
    p = plots[key]; px, py = ylab_w + (pw - p.width) // 2, header_h + (ph - p.height) // 2
    cell.paste(p, (px, py))
    yt = Image.new("RGBA", (ph, ylab_w), (255, 255, 255, 0))
    ImageDraw.Draw(yt).text((ph // 2, ylab_w // 2 + 3), yl, fill="black", font=fax, anchor="mm")
    cell.paste(yt.rotate(90, expand=True), (0, header_h), yt.rotate(90, expand=True))
    dd.text((ylab_w + pw // 2, header_h + ph + 8), xl, fill="black", font=fax, anchor="ma")
    return cell

def fit_w(im, w): return im.resize((w, round(w * im.height / im.width)))
def fit_h(im, h): return im.resize((round(h * im.width / im.height), h))

photos = Image.open(f"{FIG}/FigS_pdl1.png").convert("RGB")     # PD-L1 IHC images
ihc_q  = Image.open(f"{FIG}/Fig_quant_PDL1.png").convert("RGB")  # PD-L1 IHC quantification
pdl1_g = Image.open(f"{FIG}/Fig_pdl1_flow.png").convert("RGB")   # PD-L1 flow time course
pd1_g  = Image.open(f"{FIG}/Fig_pd1_flow.png").convert("RGB")    # PD-1 flow time course
gate_cells = {k: layout_cell(k) for k in ("PDL1", "PD1")}

W, colgap, rowgap, lm, margin = 4200, 130, 80, 150, 60
# top IHC row: photos (left) + quant (right), vertically centred to a common height
Pw = 2520; photos_r = fit_w(photos, Pw)
quant_r = fit_w(ihc_q, W - Pw - colgap)
Ht = max(photos_r.height, quant_r.height)
# flow rows: gate (left) + line graph (right); pick a height so gate + graph span W
ca = gate_cells["PDL1"].width / gate_cells["PDL1"].height
la = pdl1_g.width / pdl1_g.height
Hf = round((W - colgap) / (ca + la))

def flow_row(gkey, graph):
    return fit_h(gate_cells[gkey], Hf), fit_h(graph, Hf)

gC, lC = flow_row("PDL1", pdl1_g)
gD, lD = flow_row("PD1", pd1_g)

Wtot = lm + W + margin
Htot = margin + Ht + rowgap + Hf + rowgap + Hf + margin
canvas = Image.new("RGB", (Wtot, Htot), "white"); d = ImageDraw.Draw(canvas)
fl = font(120)

# --- top IHC row ---
yT = margin
canvas.paste(photos_r, (lm, yT + (Ht - photos_r.height) // 2))
xq = lm + Pw + colgap
canvas.paste(quant_r, (xq, yT + (Ht - quant_r.height) // 2))
d.text((30, yT), "A", fill="black", font=fl, anchor="la")
d.text((xq - 96, yT), "B", fill="black", font=fl, anchor="la")

# --- flow row C (PD-L1) ---
yC = yT + Ht + rowgap
canvas.paste(gC, (lm, yC)); canvas.paste(lC, (lm + gC.width + colgap, yC))
d.text((30, yC), "C", fill="black", font=fl, anchor="la")
# --- flow row D (PD-1) ---
yD = yC + Hf + rowgap
canvas.paste(gD, (lm, yD)); canvas.paste(lD, (lm + gD.width + colgap, yD))
d.text((30, yD), "D", fill="black", font=fl, anchor="la")

out = f"{FIG}/Fig_checkpoints_full.png"
canvas.save(out, dpi=(300, 300))
print("saved", out, canvas.size, "aspect", round(Wtot / Htot, 2))
