#!/usr/bin/env python3
"""Annotated H&E montage: same layout as FigS_spheroids_HE but each spheroid is overlaid with the
color-deconvolution analysis (nuclei outlined green, necrotic/acellular areas outlined red) — the
same method as preview_he_method. Output: output/figuras/FigS_spheroids_HE_annotated.png"""
import os
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import (gaussian_filter, label, binary_fill_holes, binary_closing,
                           binary_erosion, binary_dilation)

FT = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES"
HE = os.path.join(FT, "FOTOS H&E SCANER")
HE_MICRO = os.path.join(FT, "fotos H&E MICROSCOPIO")
FIG = "output/figuras"
CONDS = ["NI_Ctrl", "NI_Rest", "NI_Act", "INF_Ctrl", "INF_Rest", "INF_Act"]
CONDLAB = ["Control", "Resting", "Activated", "Control", "Resting", "Activated"]
TIMES = ["24", "48", "96"]

# Ruifrok-Johnston H&E stain matrix (as in he_analysis.py / preview_he_method)
M = np.array([[0.65, 0.70, 0.29], [0.07, 0.99, 0.11], [0.27, 0.57, 0.78]])
M = M / np.linalg.norm(M, axis=1, keepdims=True); Minv = np.linalg.inv(M)
ACELL_THR = 0.13

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

def otsu(v):
    h, e = np.histogram(v, bins=128); p = h / h.sum(); w = np.cumsum(p)
    m = np.cumsum(p * ((e[:-1] + e[1:]) / 2)); mt = m[-1]
    sb = (mt * w - m) ** 2 / (w * (1 - w) + 1e-9); return ((e[:-1] + e[1:]) / 2)[np.nanargmax(sb)]

def largest_cc(mask):
    lbl, n = label(mask)
    if n == 0: return mask
    s = np.bincount(lbl.ravel()); s[0] = 0; return lbl == s.argmax()

def rm_small(mask, minpx):
    lbl, n = label(mask); s = np.bincount(lbl.ravel())
    keep = np.where(s >= minpx)[0]; keep = keep[keep != 0]; return np.isin(lbl, keep)

def outline(mask, it=1):
    return binary_dilation(mask, iterations=it) & ~binary_erosion(mask, iterations=it)

def annotate(rgb):
    """rgb: HxWx3 float in [0,255] (white-balanced). Returns annotated uint8, plus (%nuclei, %necrotic)."""
    s = min(rgb.shape[:2]); f = s / 1100.0            # scale detection params to image size
    sig = max(4, 16 * f); er = max(4, round(16 * f)); cl = max(3, round(11 * f)) | 1
    a = np.clip(rgb / 255.0, 1e-6, 1.0); od = -np.log(a)
    Hc = (od.reshape(-1, 3) @ Minv).reshape(rgb.shape)[..., 0]
    stained = od.sum(2) > 0.25
    if stained.sum() < 50:
        return rgb.astype("uint8"), (0.0, 0.0)
    thr = max(otsu(Hc[stained]), Hc[stained].mean() + 0.5 * Hc[stained].std())
    nuclei0 = stained & (Hc > thr)
    dens = gaussian_filter(nuclei0.astype(float), sigma=sig)
    body = binary_fill_holes(largest_cc(binary_closing(dens > 0.05, np.ones((cl, cl)), iterations=2)))
    nuclei = body & nuclei0
    inner = binary_erosion(body, iterations=er)
    acell = rm_small(inner & (dens < ACELL_THR), round(200 * f * f))
    acell = binary_fill_holes(acell)      # count the whole low-density zone; do NOT subtract the
                                          # few nuclei that fall inside a red area
    ov = rgb.copy()
    ov[outline(acell, max(1, round(2 * f)))] = [220, 0, 0]     # low-density area, red
    ov[outline(nuclei, 1)] = [0, 140, 0]                        # nuclei green
    tot = max(body.sum(), 1)
    return np.clip(ov, 0, 255).astype("uint8"), (100 * nuclei.sum() / tot, 100 * acell.sum() / tot)


def annotate_disagg(rgb):
    """Disaggregated-spheroid mode (loss of integrity): the spheroid has fallen apart into loose,
    largely pyknotic cell clusters. Only the compact, cohesive high-density tissue counts as viable
    (green); everything else inside the spheroid footprint — loose cells, gaps and eosinophilic
    necrotic material — is low-density / necrotic (red). Used for the two INF · Activated spheroids
    (48 h, 96 h) whose integrity loss exceeds 70 %."""
    a = np.clip(rgb / 255.0, 1e-6, 1.0); od = -np.log(a)
    Hc = (od.reshape(-1, 3) @ Minv).reshape(rgb.shape)[..., 0]
    stained = od.sum(2) > 0.18
    if stained.sum() < 50:
        return rgb.astype("uint8"), (0.0, 0.0)
    thr = max(otsu(Hc[stained]), Hc[stained].mean() + 0.5 * Hc[stained].std())
    nuclei0 = stained & (Hc > thr)
    # body = the coherent spheroid, from the smoothed tissue-density field: sparse scattered debris
    # smooths below threshold and is dropped, so the selection never reaches out to the background.
    tf = gaussian_filter(stained.astype(float), sigma=16)
    body = binary_fill_holes(largest_cc(tf > 0.55))
    tissue = body & stained                                 # actual tissue inside the body (no background)
    dens = gaussian_filter(nuclei0.astype(float), sigma=8)
    viable = binary_fill_holes(binary_closing(rm_small(dens > 0.55, 500), np.ones((9, 9)))) & tissue
    nec = tissue & ~viable
    # red is painted only on tissue pixels (never on background/empty gaps); viable tissue outlined green
    ov = rgb.astype(float); red = np.zeros_like(ov); red[..., 0] = 220.0
    ov = np.where(nec[..., None], 0.6 * ov + 0.4 * red, ov)
    ov = ov.astype("uint8")
    ov[outline(viable, 2)] = [0, 140, 0]
    tot = max(tissue.sum(), 1)
    return np.clip(ov, 0, 255).astype("uint8"), (100 * (tissue & viable).sum() / tot, 100 * nec.sum() / tot)

def fit_square(path, box):
    im = Image.open(path).convert("RGB"); w, h = im.size; s = min(w, h)
    l, t = (w - s) // 2, (h - s) // 2
    return np.asarray(im.crop((l, t, l + s, t + s)).resize((box, box))).astype(float)

def build():
    idx = index(HE, (".png",))
    fb = {k: v for k, v in index(HE_MICRO, (".tif",)).items() if "20x" in v.lower()}
    lg, cell, gap = 210, 620, 14
    envh, condh, margin = 62, 56, 40
    grid_w = 6 * cell + 5 * gap
    W = margin + lg + grid_w + margin
    y_env = margin + 70; y_cond = y_env + envh; y0 = y_cond + condh
    H = y0 + 3 * cell + 2 * gap + margin + 60
    canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
    def colx(i): return margin + lg + i * (cell + gap)
    d.text((W // 2, margin), "Spheroid histology (H&E) — nuclei (green) / necrotic areas (red)",
           fill="black", font=font(50), anchor="ma")
    d.text((colx(0) + (3 * cell + 2 * gap) // 2, y_env), "Basal (non-inflammatory)", fill="black", font=font(48), anchor="ma")
    d.text((colx(3) + (3 * cell + 2 * gap) // 2, y_env), "Inflammatory", fill="black", font=font(48), anchor="ma")
    for i, c in enumerate(CONDLAB):
        d.text((colx(i) + cell // 2, y_cond + 2), c, fill="black", font=font(40), anchor="ma")
    xdiv = colx(3) - gap // 2
    d.line([(xdiv, y_cond), (xdiv, y0 + 3 * cell + 2 * gap)], fill=(150, 150, 150), width=3)
    rows = []
    envmap = {"NI": "Basal", "INF": "Inflammatory"}
    cmap = {"Ctrl": "Control", "Rest": "Resting", "Act": "Activated"}
    for r, t in enumerate(TIMES):
        y = y0 + r * (cell + gap)
        d.text((margin + lg - 16, y + cell // 2), f"{t} h", fill="black", font=font(48), anchor="rm")
        for i, c in enumerate(CONDS):
            x = colx(i); f = idx.get((t, c)); src = HE
            if not f and (t, c) in fb: f, src = fb[(t, c)], HE_MICRO
            if f:
                arr = white_balance(fit_square(os.path.join(src, f), cell))
                if (t, c) in (("48", "INF_Act"), ("96", "INF_Act")):   # disaggregated spheroids
                    ov, (pn, pa) = annotate_disagg(arr)
                else:
                    ov, (pn, pa) = annotate(arr)
                env, cond = envmap[c.split("_")[0]], cmap[c.split("_")[1]]
                rows.append((int(t), env, cond, round(pa, 2)))
                canvas.paste(Image.fromarray(ov), (x, y))
                d.rectangle([x, y, x + cell - 1, y + cell - 1], outline=(0, 0, 0), width=2)
            else:
                d.rectangle([x, y, x + cell, y + cell], fill=(238, 238, 238), outline=(200, 200, 200), width=2)
                d.text((x + cell // 2, y + cell // 2), "n/a", fill=(140, 140, 140), font=font(46), anchor="mm")
    # legend
    ly = y0 + 3 * cell + 2 * gap + 18; lx = margin + lg
    d.rectangle([lx, ly + 6, lx + 42, ly + 30], outline=(0, 140, 0), width=6)
    d.text((lx + 56, ly + 2), "nuclei", fill="black", font=font(38))
    d.rectangle([lx + 300, ly + 6, lx + 342, ly + 30], outline=(220, 0, 0), width=6)
    d.text((lx + 356, ly + 2), "necrotic / acellular area", fill="black", font=font(38))
    canvas.save(os.path.join(FIG, "FigS_spheroids_HE_annotated.png"), dpi=(300, 300))
    print("saved FigS_spheroids_HE_annotated.png", canvas.size)
    os.makedirs("output/tidy", exist_ok=True)
    with open("output/tidy/he_lowdensity.csv", "w") as fh:
        fh.write("time,environment,condition,pct_low_density\n")
        for tt, env, cond, pa in sorted(rows, key=lambda z: (z[1], z[2], z[0])):
            fh.write(f"{tt},{env},{cond},{pa}\n")
    print("saved output/tidy/he_lowdensity.csv")

build()
