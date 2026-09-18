#!/usr/bin/env python3
"""Ki-67 figure, unified layout: left = 2 x 3 IHC image grid (Basal|Inflammatory x
Control/Resting/Activated), right-top = Ki-67 quantification chart.
Output: output/figuras/Fig_IHC_Ki67.png"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from PIL import Image, ImageDraw, ImageFont
from img_norm import fit_object
IH = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES/FOTOS IH"
FIG = "output/figuras"
# grid order: rows = Control/Resting/Activated ; cols = Basal, Inflammatory
IMGS = {
    ("Control", "Basal"):        "1. Esf no act no inf 48 hrs Ki67.png",
    ("Resting", "Basal"):        "2. Esf +PBMC no act no inf 48 hrs Ki67.png",
    ("Activated", "Basal"):      "3. Esf +PBMC act no inf 48 hrs Ki67.png",
    ("Control", "Inflammatory"): "4.Esf inf 48 hrs Ki67.png",
    ("Resting", "Inflammatory"): "5. Esf + PBMC no act inf 48 hrs Ki67.png",
    ("Activated", "Inflammatory"): "6. Esf + PBMC act inf 24 hrs Ki67.png",
}
ROWS = ["Control", "Resting", "Activated"]
COLS = ["Basal", "Inflammatory"]
def font(sz, bold=True):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans%s.ttf" % ("-Bold" if bold else ""), sz)

cell, gap, tgut, condh, margin, colgap = 680, 16, 250, 60, 40, 80
x0 = margin + 60
grid_w = 2 * cell + gap
grid_h = 3 * cell + 2 * gap
y_cond = margin; y_grid = y_cond + condh
def colx(i): return x0 + tgut + i * (cell + gap)

# graph top-right
G = Image.open(f"{FIG}/Fig_quant_Ki67.png").convert("RGB")
Gw = 1420; Gh = round(Gw * G.height / G.width); G = G.resize((Gw, Gh))
xG = x0 + tgut + grid_w + colgap
W = xG + Gw + margin
H = margin + condh + grid_h + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)

fcond, frow, flet = font(46), font(46), font(64)
for i, c in enumerate(COLS):
    d.text((colx(i) + cell // 2, y_cond + 2), c, fill="black", font=fcond, anchor="ma")
for r, rlab in enumerate(ROWS):
    y = y_grid + r * (cell + gap)
    d.text((x0 + tgut - 18, y + cell // 2), rlab, fill="black", font=frow, anchor="rm")
    for i, cenv in enumerate(COLS):
        p = Image.open(os.path.join(IH, IMGS[(rlab, cenv)])).convert("RGB")
        canvas.paste(fit_object(p, cell, cell), (colx(i), y))

# graph: top-right, top-aligned with the grid
d.text((10, y_grid - 4), "A", fill="black", font=flet, anchor="lb")
d.text((xG - 54, y_grid - 4), "B", fill="black", font=flet, anchor="lb")
canvas.paste(G, (xG, y_grid))

canvas.save(f"{FIG}/Fig_IHC_Ki67.png", dpi=(300, 300))
print("saved Fig_IHC_Ki67.png", canvas.size, "aspect", round(W / H, 2))
