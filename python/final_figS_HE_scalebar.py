#!/usr/bin/env python3
"""FigS_spheroids_HE: each photo is CONTAINED whole inside its square cell (no square-crop), so the
spheroid keeps the same framing as the source file and never spills past the border. Empty area is
padded white (matches the H&E background). A 100 µm scale bar is burned into each cell, black and
without a number in the Fig1A style (8 px thick, 30 px margin, bottom-right). All bars are identical
in length (133 px = the 100 µm reference, same proportion as Fig_HE) and sit at the same height in
every cell. Output: output/figuras/FigS_spheroids_HE.png"""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import (gaussian_filter, label, binary_fill_holes, binary_closing,
                           binary_dilation, binary_erosion)

FT = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES"
HE = os.path.join(FT, "FOTOS H&E SCANER")
HE_MICRO = os.path.join(FT, "fotos H&E MICROSCOPIO")
FIG = "output/figuras"
CONDS = ["NI_Ctrl", "NI_Rest", "NI_Act", "INF_Ctrl", "INF_Rest", "INF_Act"]
CONDLAB = ["Control", "Resting", "Activated", "Control", "Resting", "Activated"]
TIMES = ["24", "48", "96"]
CELL, GAP = 620, 14
BAR_TH, BM = 8, 30                   # Fig1A bar style (thickness / margin from content edges)
GRID_W = 2                           # thin black border drawn around every image cell (grid)
BAR_PX = 133                         # 100 µm reference, uniform in every cell (Fig_HE proportion)

# Ruifrok-Johnston H&E stain matrix (hematoxylin channel = nuclei), as in compose_HE_annotated.py
M = np.array([[0.65, 0.70, 0.29], [0.07, 0.99, 0.11], [0.27, 0.57, 0.78]])
M = M / np.linalg.norm(M, axis=1, keepdims=True); Minv = np.linalg.inv(M)


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
    return {k: sorted(v, key=len)[0] for k, v in g.items()}


def white_balance(a, target=245.0):
    ref = np.clip(np.percentile(a.reshape(-1, 3), 92, axis=0), 60, None)
    return np.clip(a * np.clip(target / ref, 0.7, 1.6), 0, 255)


def largest_cc(mask):
    lbl, n = label(mask)
    if n == 0:
        return mask
    s = np.bincount(lbl.ravel()); s[0] = 0
    return lbl == s.argmax()


def otsu(v):
    h, e = np.histogram(v, bins=128); p = h / h.sum(); w = np.cumsum(p)
    m = np.cumsum(p * ((e[:-1] + e[1:]) / 2)); mt = m[-1]
    sb = (mt * w - m) ** 2 / (w * (1 - w) + 1e-9)
    return ((e[:-1] + e[1:]) / 2)[np.nanargmax(sb)]


def whiten_bg(rgb):
    """Keep only the nucleated spheroid body and set everything around it (scanner background, tissue
    debris, acellular pink/eosin folds) to pure white. The body is the largest blob of high NUCLEAR
    (hematoxylin) density, so eosin-only folds adhering to the spheroid are excluded even when they
    touch it. Same nuclear-density body used by the annotated H&E montage."""
    hh, ww = rgb.shape[:2]; s = min(hh, ww); f = s / 1100.0
    sig = max(4, 16 * f); cl = max(3, round(11 * f)) | 1
    a = np.clip(rgb / 255.0, 1e-6, 1.0); od = -np.log(a)
    Hc = (od.reshape(-1, 3) @ Minv).reshape(rgb.shape)[..., 0]
    stained = od.sum(2) > 0.25
    if stained.sum() < 50:
        return rgb.astype("uint8")
    thr = max(otsu(Hc[stained]), Hc[stained].mean() + 0.5 * Hc[stained].std())
    nuclei0 = stained & (Hc > thr)
    dens = gaussian_filter(nuclei0.astype(float), sigma=sig)
    body = binary_fill_holes(largest_cc(binary_closing(dens > 0.05, np.ones((cl, cl)), iterations=2)))
    if body.sum() < 50:
        return rgb.astype("uint8")
    # opening: erode -> largest blob -> dilate back. Severs thin full-width tissue bands that touch
    # the spheroid (e.g. INF Control 48 h top strip) while the thick spheroid core survives.
    op = 20                                                          # content's long side is always 620
    core = largest_cc(binary_erosion(body, iterations=op))
    if core.sum() >= 50:
        body = binary_fill_holes(binary_dilation(core, iterations=op))
    body = binary_dilation(body, iterations=max(2, round(6 * f)))    # include cytoplasmic rim
    out = rgb.copy(); out[~body] = 255
    return out.astype("uint8")


def fit(path, box):
    """Contain the whole image inside a white box square, then draw the 100 µm bar (scaled by how
    much this photo was shrunk) at the bottom-right of the actual photo content."""
    im = Image.open(path).convert("RGB"); w, h = im.size
    scale = box / max(w, h)
    nw, nh = max(1, round(w * scale)), max(1, round(h * scale))
    content = Image.fromarray(whiten_bg(white_balance(np.asarray(im.resize((nw, nh))).astype(float))))
    cell = Image.new("RGB", (box, box), "white")
    ox, oy = (box - nw) // 2, (box - nh) // 2
    cell.paste(content, (ox, oy))
    d = ImageDraw.Draw(cell)
    xr = box - BM; xl = xr - BAR_PX; yb = box - BM    # uniform length, fixed corner -> identical bars
    d.rectangle([xl, yb - BAR_TH, xr, yb], fill="black")
    return cell


idx = index(HE, (".png",))
fb = {k: v for k, v in index(HE_MICRO, (".tif",)).items() if "20x" in v.lower()}
lg, cell, gap = 210, CELL, GAP
envh, condh, margin = 62, 56, 40
grid_w = 6 * cell + 5 * gap
W = margin + lg + grid_w + margin
y_env = margin + 70; y_cond = y_env + envh; y0 = y_cond + condh
H = y0 + 3 * cell + 2 * gap + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)


def colx(i): return margin + lg + i * (cell + gap)


d.text((W // 2, margin), "Spheroid histology (H&E)", fill="black", font=font(58), anchor="ma")
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
        x = colx(i); f = idx.get((t, c)); src = HE
        if not f and (t, c) in fb:
            f, src = fb[(t, c)], HE_MICRO
        if f:
            canvas.paste(fit(os.path.join(src, f), cell), (x, y))
        else:
            d.rectangle([x, y, x + cell, y + cell], fill=(255, 255, 255))
            d.text((x + cell // 2, y + cell // 2), "n/a", fill=(140, 140, 140), font=font(46), anchor="mm")
        d.rectangle([x, y, x + cell - 1, y + cell - 1], outline="black", width=GRID_W)   # thin grid
canvas.save(os.path.join(FIG, "FigS_spheroids_HE.png"), dpi=(300, 300))
print("saved FigS_spheroids_HE.png", canvas.size)
