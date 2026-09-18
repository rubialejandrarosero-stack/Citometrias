#!/usr/bin/env python3
"""Supplementary spheroid figures: all timepoints (24/48/96 h) x 6 conditions
(Basal|Inflammatory x Control/Resting/Activated), one figure for brightfield and one for H&E.
Uniform square cells. Outputs: FigS_spheroids_brightfield.png / FigS_spheroids_HE.png"""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from img_norm import fit_object
FT = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES"
BF = os.path.join(FT, "FOTOS ESFEROIDES CONTROL"); HE = os.path.join(FT, "FOTOS H&E SCANER")
HE_MICRO = os.path.join(FT, "fotos H&E MICROSCOPIO")   # fallback for cells missing in the scanner set
FIG = "output/figuras"
CONDS = ["NI_Ctrl", "NI_Rest", "NI_Act", "INF_Ctrl", "INF_Rest", "INF_Act"]
CONDLAB = ["Control", "Resting", "Activated", "Control", "Resting", "Activated"]
TIMES = ["24", "48", "96"]

def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

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
    return {k: sorted(v, key=len)[0] for k, v in g.items()}   # pick the simplest name

def fit(path, box):
    im = Image.open(path).convert("RGB"); w, h = im.size; s = min(w, h)
    l, t = (w - s) // 2, (h - s) // 2
    return im.crop((l, t, l + s, t + s)).resize((box, box))

def white_balance(im, target=245.0):
    """Unify the background: map each channel's bright point (the background) to near-white,
    so all H&E backgrounds look consistent. Uses a per-channel high percentile (robust)."""
    a = np.asarray(im).astype(float)
    ref = np.clip(np.percentile(a.reshape(-1, 3), 92, axis=0), 60, None)   # background colour
    scale = np.clip(target / ref, 0.7, 1.6)
    return Image.fromarray(np.clip(a * scale, 0, 255).astype("uint8"))

def build(folder, exts, out, title, fallback=None, wb=False):
    idx = index(folder, exts)
    fb = index(fallback, (".tif",)) if fallback else {}
    fb = {k: v for k, v in fb.items() if "20x" in v.lower()}   # prefer the 20X overview
    lg, cell, gap = 210, 620, 14
    envh, condh, margin, rowlab = 62, 56, 40, 210
    grid_w = 6 * cell + 5 * gap
    W = margin + lg + grid_w + margin
    y_env = margin + 70
    y_cond = y_env + envh
    y0 = y_cond + condh
    H = y0 + 3 * cell + 2 * gap + margin
    canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
    def colx(i): return margin + lg + i * (cell + gap)
    d.text((W // 2, margin), title, fill="black", font=font(58), anchor="ma")
    d.text((colx(0) + (3 * cell + 2 * gap) // 2, y_env), "Basal (non-inflammatory)", fill="black", font=font(48), anchor="ma")
    d.text((colx(3) + (3 * cell + 2 * gap) // 2, y_env), "Inflammatory", fill="black", font=font(48), anchor="ma")
    for i, c in enumerate(CONDLAB):
        d.text((colx(i) + cell // 2, y_cond + 2), c, fill="black", font=font(40), anchor="ma")
    xdiv = colx(3) - gap // 2
    d.line([(xdiv, y_cond), (xdiv, y0 + 3 * cell + 2 * gap)], fill=(150, 150, 150), width=3)
    for r, t in enumerate(TIMES):
        y = y0 + r * (cell + gap)
        d.text((margin + lg - 16, y + cell // 2), f"{t} h", fill="black", font=font(48), anchor="rm")
        for i, c in enumerate(CONDS):
            x = colx(i)
            f = idx.get((t, c))
            src = folder
            if not f and (t, c) in fb:                # borrow missing cell from the microscope set
                f, src = fb[(t, c)], fallback
            if f:
                img = fit(os.path.join(src, f), cell)
                if wb:                       # H&E: unify background to white (no size change)
                    img = white_balance(img)
                canvas.paste(img, (x, y))
            else:
                d.rectangle([x, y, x + cell, y + cell], fill=(238, 238, 238), outline=(200, 200, 200), width=2)
                d.text((x + cell // 2, y + cell // 2), "n/a", fill=(140, 140, 140), font=font(46), anchor="mm")
    canvas.save(os.path.join(FIG, out), dpi=(300, 300))
    print("saved", out, canvas.size)

build(BF, (".jpg", ".jpeg"), "FigS_spheroids_brightfield.png", "Spheroid morphology (brightfield)")
build(HE, (".png",), "FigS_spheroids_HE.png", "Spheroid histology (H&E)", fallback=HE_MICRO, wb=True)
