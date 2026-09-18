#!/usr/bin/env python3
"""Figure 1A (final): brightfield spheroid montage, 6 conditions x 3 times, NO titles/labels.
A uniform scale bar (identical to the Basal Control 24 h one: 73 px long, 8 px thick) is drawn in
the bottom-right corner of every cell at the same margin; the burned-in microscope bar is covered.
Output: output/figuras_finales_11082026/Fig1A_brightfield.png"""
import os
import numpy as np
from PIL import Image, ImageDraw

BF = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES/FOTOS ESFEROIDES CONTROL"
OUT = "output/figuras_finales_11082026"; os.makedirs(OUT, exist_ok=True)
CONDS = ["NI_Ctrl", "NI_Rest", "NI_Act", "INF_Ctrl", "INF_Rest", "INF_Act"]
TIMES = ["24", "48", "96"]
CELL, GAP, MARGIN = 620, 14, 20
BAR_PX, BAR_TH, BM = 73, 8, 30      # bar length / thickness / margin from right & bottom edges

def classify(name):
    n = name.lower()
    t = "24" if "24" in n else "48" if "48" in n else "96" if "96" in n else "?"
    inf = "NI" if ("no inf" in n or "noinf" in n or "no_inf" in n) else "INF"
    if "ctl" in n or "control" in n or ("pbmc" not in n and "esf" in n):
        cond = "Ctrl"
    elif "no act" in n or "noact" in n or "no_act" in n:
        cond = "Rest"
    elif "act" in n:
        cond = "Act"
    else:
        cond = "Ctrl"
    return t, f"{inf}_{cond}"

def index(folder, exts):
    g = {}
    for f in sorted(os.listdir(folder)):
        if f.lower().endswith(exts) and "original" not in f.lower():
            g.setdefault(classify(f), []).append(f)
    return {k: sorted(v, key=len)[0] for k, v in g.items()}

def fit(path, box):
    im = Image.open(path).convert("RGB"); w, h = im.size; s = min(w, h)
    l, t = (w - s) // 2, (h - s) // 2
    return im.crop((l, t, l + s, t + s)).resize((box, box))

def uniform_bar(img):
    """Detect and erase the burned-in bottom-right bar (fill with local background), then draw the
    uniform bar (same length as the Basal Control 24 h reference)."""
    a = np.array(img); g = a.mean(2)
    ry0, rx0 = int(CELL * 0.80), int(CELL * 0.42)      # bottom-right search region
    dark = g < 75
    best = None
    for y in range(ry0, CELL):
        xs = np.where(dark[y, rx0:])[0]
        if len(xs) < 25:
            continue
        segs = np.split(xs, np.where(np.diff(xs) > 3)[0] + 1)
        seg = max(segs, key=len); run = seg.max() - seg.min()
        if run > 40 and (best is None or run > best[0]):
            best = (run, y, seg.min() + rx0, seg.max() + rx0)
    if best:
        run, y, x0, x1 = best
        xc = (x0 + x1) // 2; ys = np.where(dark[ry0:, xc])[0] + ry0
        by0, by1 = ys.min(), ys.max(); pad = 5
        bg = np.median(a[max(0, by0 - 26):by0 - 5, x0:x1].reshape(-1, 3), axis=0).astype("uint8")
        a[max(0, by0 - pad):by1 + pad, max(0, x0 - pad):min(CELL, x1 + pad)] = bg
        img = Image.fromarray(a)
    d = ImageDraw.Draw(img)
    xr = CELL - BM; xl = xr - BAR_PX; yb = CELL - BM
    d.rectangle([xl, yb - BAR_TH, xr, yb], fill="black")
    return img

idx = index(BF, (".jpg", ".jpeg"))
grid_w = 6 * CELL + 5 * GAP
grid_h = 3 * CELL + 2 * GAP
W = MARGIN * 2 + grid_w
H = MARGIN * 2 + grid_h
canvas = Image.new("RGB", (W, H), "white")
def colx(i): return MARGIN + i * (CELL + GAP)
for r, t in enumerate(TIMES):
    y = MARGIN + r * (CELL + GAP)
    for i, c in enumerate(CONDS):
        f = idx.get((t, c))
        if f:
            canvas.paste(uniform_bar(fit(os.path.join(BF, f), CELL)), (colx(i), y))
        else:
            ImageDraw.Draw(canvas).rectangle([colx(i), y, colx(i) + CELL, y + CELL],
                                             fill=(238, 238, 238), outline=(200, 200, 200), width=2)
out = os.path.join(OUT, "Fig1A_brightfield.png")
canvas.save(out, dpi=(300, 300))
print("saved", out, canvas.size)
