#!/usr/bin/env python3
"""Unified Figure 1 (morphology): row A brightfield + row B H&E of representative 96 h
spheroids across the 6 conditions (Basal|Inflammatory x Control/Resting/Activated), then
panels C/D/E = area / diameter / circularity plots. Uses the user's processed images.
Output: output/figuras/Fig1_morphology_full.png"""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from img_norm import fit_object
FT = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES"
BF = os.path.join(FT, "FOTOS ESFEROIDES CONTROL")
HE = os.path.join(FT, "FOTOS H&E SCANER")
FIG = "output/figuras"
CONDS = ["Control", "Resting", "Activated", "Control", "Resting", "Activated"]
BF_F = ["96 NO INF-CTL.jpg", "96 NO INF NO ACT.jpg", "96 NO INF ACT.jpg",
        "96 INF CTL.jpg", "96 INF NO ACT.jpg", "96 INF ACT.jpg"]
HE_F = ["Esf no inf 96 hrs.png", "Esf+PBMC no act no inf 96 hrs.png", "Esf+PBMC Act no inf 96 hrs.png",
        "Esf inf 96 hrs.png", "Esf+ PBMC no act inf 96 hrs.png", "Esf+PBMC Act Inf 96 hrs.png"]
PLOTS = ["Fig1d_area.png", "Fig1d_diameter.png", "Fig1d_circularity.png"]

def font(sz, bold=True):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-%s.ttf" % ("Bold" if bold else "Oblique"), sz)

def fit(path, box):
    """Center-crop to a square and resize, so every cell is the same size and fully filled."""
    im = Image.open(path).convert("RGB")
    w, h = im.size; s = min(w, h)
    l, t = (w - s) // 2, (h - s) // 2
    return im.crop((l, t, l + s, t + s)).resize((box, box))

def white_balance(im, target=245.0):
    """Unify the background to near-white (per-channel bright reference). Used for H&E only."""
    a = np.asarray(im).astype(float)
    ref = np.clip(np.percentile(a.reshape(-1, 3), 92, axis=0), 60, None)
    scale = np.clip(target / ref, 0.7, 1.6)
    return Image.fromarray(np.clip(a * scale, 0, 255).astype("uint8"))

lg, cell, gap = 250, 760, 16                       # left gutter, image cell, gap
envh, condh = 66, 60                               # header strip heights
grid_w = 6 * cell + 5 * gap
W = lg + grid_w
plot_w = (grid_w - 2 * gap) // 3
plot_h = round(plot_w * 1980 / 5100)
margin, rowgap, plotgap = 40, 16, 170
y_env = margin
y_cond = y_env + envh
y_A = y_cond + condh
y_B = y_A + cell + rowgap
y_plots = y_B + cell + plotgap
H = y_plots + plot_h + margin
canvas = Image.new("RGB", (W + margin, H), "white"); d = ImageDraw.Draw(canvas)

def colx(i): return lg + i * (cell + gap)

# environment spanning labels + divider
fenv, fcond, frow, flet = font(56), font(46), font(50), font(64)
d.text((colx(0) + (3 * cell + 2 * gap) // 2, y_env), "Basal (non-inflammatory)", fill="black", font=fenv, anchor="ma")
d.text((colx(3) + (3 * cell + 2 * gap) // 2, y_env), "Inflammatory", fill="black", font=fenv, anchor="ma")
xdiv = colx(3) - gap // 2
d.line([(xdiv, y_cond), (xdiv, y_B + cell)], fill=(150, 150, 150), width=3)
for i, c in enumerate(CONDS):
    d.text((colx(i) + cell // 2, y_cond + 2), c, fill="black", font=fcond, anchor="ma")

# row A brightfield, row B H&E
for i, f in enumerate(BF_F):
    canvas.paste(fit(os.path.join(BF, f), cell), (colx(i), y_A))
for i, f in enumerate(HE_F):
    canvas.paste(white_balance(fit(os.path.join(HE, f), cell)), (colx(i), y_B))
d.text((lg - 16, y_A + cell // 2), "Brightfield", fill="black", font=frow, anchor="rm")
d.text((lg - 16, y_B + cell // 2), "H&E", fill="black", font=frow, anchor="rm")

# panels C/D/E
for i, p in enumerate(PLOTS):
    im = Image.open(os.path.join(FIG, p)).convert("RGB").resize((plot_w, plot_h))
    x = lg + i * (plot_w + gap)
    canvas.paste(im, (x, y_plots))
    d.text((x + 4, y_plots - 45), "CDE"[i], fill="black", font=flet, anchor="lb")  # in the gap above each plot

d.text((10, y_A - 4), "A", fill="black", font=flet, anchor="lb")
d.text((10, y_B - 4), "B", fill="black", font=flet, anchor="lb")
canvas.save(os.path.join(FIG, "Fig1_morphology_full.png"), dpi=(300, 300))
print("saved Fig1_morphology_full.png", canvas.size)
