#!/usr/bin/env python3
"""Dedicated H&E histology figure: representative sections, Control vs Activated (inflammatory)
at 24/48/96 h. (Necrosis-annotated panel to be added.) Output: output/figuras/Fig_HE.png"""
import os, numpy as np
from PIL import Image, ImageDraw, ImageFont
HE = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES/FOTOS H&E SCANER"
FIG = "output/figuras"
ROWS = [
    ("24 h", "Esf inf 24 hrs.png", "Esf +PBMC act inf 24 hrs.png"),
    ("48 h", "Esf inf 48 hrs.png", "Esf+PBMC Act Inf 48 hrs.png"),
    ("96 h", "Esf inf 96 hrs.png", "Esf+PBMC Act inf 96 hrs.png"),
]
def font(sz, bold=True):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans%s.ttf" % ("-Bold" if bold else ""), sz)

def white_balance(im, target=245.0):
    a = np.asarray(im).astype(float)
    ref = np.clip(np.percentile(a.reshape(-1, 3), 92, axis=0), 60, None)
    return Image.fromarray(np.clip(a * np.clip(target / ref, 0.7, 1.6), 0, 255).astype("uint8"))

def fit(path, box):
    im = Image.open(path).convert("RGB"); w, h = im.size; s = min(w, h)
    l, t = (w - s) // 2, (h - s) // 2
    return white_balance(im.crop((l, t, l + s, t + s)).resize((box, box)))

cell, gap, rowgap, tgut, envh, condh, margin = 700, 16, 16, 150, 66, 58, 40
x0 = margin + 60
grid_w = 2 * cell + gap
strips_h = envh + condh + 3 * cell + 2 * rowgap
# panel B (annotated image) on the right
Bimg = Image.open(f"{FIG}/he_both_act24.png").convert("RGB")
Bw = 1320; Bh = round(Bw * Bimg.height / Bimg.width); Bimg = Bimg.resize((Bw, Bh))
colgap = 150                              # wider gutter so panel B and its letter are not cramped
xB = x0 + tgut + grid_w + colgap
W = xB + Bw + margin
H = margin + strips_h + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)

y_env = margin; y_cond = y_env + envh; y_grid = y_cond + condh
def colx(i): return x0 + tgut + i * (cell + gap)
fenv, fcond, ftime, flet, fscale = font(48), font(44), font(40), font(64), font(34, False)
d.text((x0 + tgut + grid_w // 2, y_env), "Inflammatory", fill="black", font=fenv, anchor="ma")
d.text((colx(0) + cell // 2, y_cond + 2), "Control", fill="black", font=fcond, anchor="ma")
d.text((colx(1) + cell // 2, y_cond + 2), "Activated", fill="black", font=fcond, anchor="ma")
for r, (tlab, cf, af) in enumerate(ROWS):
    y = y_grid + r * (cell + rowgap)
    canvas.paste(fit(os.path.join(HE, cf), cell), (colx(0), y))
    canvas.paste(fit(os.path.join(HE, af), cell), (colx(1), y))
    d.text((x0 + tgut - 18, y + cell // 2), tlab, fill="black", font=ftime, anchor="rm")
# scale bar (estimated ~0.6 µm/px on the down-sized 40x scan; ~150 px = 100 µm — approximate)
BAR_PX = 150
xb = colx(1) + cell - 26; yb = y_grid + 2 * (cell + rowgap) + cell - 30
d.rectangle([xb - BAR_PX, yb - 6, xb, yb + 8], fill="white", outline="black", width=3)
d.text(((xb - BAR_PX + xb) // 2, yb - 14), "100 µm", fill="black", font=fscale, anchor="mb",
       stroke_width=3, stroke_fill="white")
d.text((10, y_grid - 4), "A", fill="black", font=flet, anchor="lb")

# ---- panel B: annotated image (nuclei + necrotic), vertically centred on the grid ----
fcap, fleg = font(40), font(34, False)
cap_h, gap_ci, gap_il = 56, 12, 30
grid_h = 3 * cell + 2 * rowgap
groupH = cap_h + gap_ci + Bh + gap_il + 40
yTop = y_grid + (grid_h - groupH) // 2
d.text((xB, yTop), "Inflammatory · Activated · 24 h", fill="black", font=fcap, anchor="la")
yImg = yTop + cap_h + gap_ci
canvas.paste(Bimg, (xB, yImg))
d.text((xB - 118, yTop + cap_h), "B", fill="black", font=flet, anchor="lb")   # clear gap from caption/image
ly = yImg + Bh + gap_il; lx = xB + 8
d.rectangle([lx, ly + 6, lx + 46, ly + 30], outline=(0, 140, 0), width=6)
d.text((lx + 60, ly + 4), "nuclei", fill="black", font=fleg, anchor="la")
lx2 = lx + 340
d.rectangle([lx2, ly + 6, lx2 + 46, ly + 30], outline=(220, 0, 0), width=6)
d.text((lx2 + 60, ly + 4), "necrotic / acellular area", fill="black", font=fleg, anchor="la")

canvas.save(f"{FIG}/Fig_HE.png", dpi=(300, 300))
print("saved Fig_HE.png", canvas.size, "aspect", round(W / H, 2))
