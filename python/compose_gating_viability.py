#!/usr/bin/env python3
"""Manuscript gating-strategy tree for the death / immune-blockade experiment (Experiment_001).
Layout (user-specified): A Live (Zombie Red) -> B CD45 split, then two rows branch off:
  top row  = PD-L1+ and EpCAM+  (both from CD45- spheroid)
  bottom   = PD-1+              (from CD3+ within CD45+)
Same treatment as FigS_gating_strategy: crop FlowJo's black export border and axis frame,
composite over white, clean axis labels, panel letters, hierarchy headers, directional arrows.
Output: output/figuras/FigS_gating_viability.png"""
import os, math, numpy as np
from PIL import Image, ImageDraw, ImageFont

LAY = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/LAYOUTS VIABILIDAD MANUSCRITO"
OUT = "output/figuras"
FDIR = "/usr/share/fonts/truetype/dejavu"
def font(sz, bold=True):
    return ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"), sz)

# key: (file, letter, title, parent, x-axis, y-axis)
panels = {
    "A": ("VIVAS",  "A", "Live cells",    "all cells",       "FSC-A", "Zombie Red"),
    "B": ("CD45-",  "B", "CD45 split",    "live cells",      "FSC-A", "CD45"),
    "C": ("PD-L1",  "C", "PD-L1+ cells",  "CD45- spheroid",  "FSC-A", "PD-L1"),
    "D": ("EpCAM",  "D", "EpCAM+ tumor",  "CD45- spheroid",  "FSC-A", "EpCAM"),
    "E": ("PD-1+ ", "E", "PD-1+ T cells", "CD3+ (CD45+)",    "CD3",   "PD-1"),
}
BLUE, ORANGE, GREY = "#3b6ea5", "#e6550d", "#555555"

def load(name):
    im = Image.open(os.path.join(LAY, name + ".png")).convert("RGBA")
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    im = Image.alpha_composite(bg, im).convert("RGB")
    g = np.array(im.convert("L"))
    r = np.where(g.max(axis=1) > 150)[0]         # trim the solid black export border
    c = np.where(g.max(axis=0) > 150)[0]
    im = im.crop((c.min(), r.min(), c.max() + 1, r.max() + 1)); g = g[r.min():r.max()+1, c.min():c.max()+1]
    H, Wd = g.shape                              # crop to the axis frame (drop FlowJo ticks + fluorophore titles)
    vc = [i for i in range(Wd) if (g[:, i] < 80).sum() > 0.6 * H]
    hr = [i for i in range(H) if (g[i, :] < 80).sum() > 0.6 * Wd]
    return im.crop((min(vc), min(hr), max(vc) + 1, max(hr) + 1))

imgs = {k: load(v[0]) for k, v in panels.items()}
pw = max(i.width for i in imgs.values()); ph = max(i.height for i in imgs.values())

ylab_w, header_h, xlab_h = 54, 82, 48
margin, title_h, gx, gy = 44, 100, 50, 58
cell_w, cell_h = ylab_w + pw, header_h + ph + xlab_h

top0 = margin + title_h
row_w = 4 * cell_w + 3 * gx                       # top row: A B C D at the same height
yTop, yBot = top0, top0 + cell_h + gy
pos = {k: (margin + i * (cell_w + gx), yTop) for i, k in enumerate(["A", "B", "C", "D"])}
pos["E"] = (margin, yBot)                          # PD-1 below Zombie Red (column A)

W = margin * 2 + row_w
Htot = yBot + cell_h + margin
canvas = Image.new("RGB", (W, Htot), "white"); d = ImageDraw.Draw(canvas)

ftitle, fhead, fpar, flet, fax, ftag = font(46), font(32), font(25, False), font(40), font(30), font(28)
d.text((W // 2, margin + 22), "Gating strategy — viability & immune checkpoint", fill="black", font=ftitle, anchor="mm")
d.text((W // 2, margin + 66), "death / immune-blockade experiment, representative sample",
       fill=(90, 90, 90), font=fpar, anchor="mm")

def draw_panel(k):
    name, letter, title, parent, xl, yl = panels[k]; cx, cy = pos[k]
    d.text((cx, cy), letter, fill="black", font=flet, anchor="la")
    d.text((cx + 52, cy + 2), title, fill="black", font=fhead, anchor="la")
    d.text((cx + 52, cy + 44), "from " + parent, fill=(90, 90, 90), font=fpar, anchor="la")
    px, py = cx + ylab_w, cy + header_h
    canvas.paste(imgs[k], (px, py))
    yt = Image.new("RGBA", (ph, ylab_w), (255, 255, 255, 0))
    ImageDraw.Draw(yt).text((ph // 2, ylab_w // 2 + 4), yl, fill="black", font=fax, anchor="mm")
    canvas.paste(yt.rotate(90, expand=True), (cx, py), yt.rotate(90, expand=True))
    d.text((px + pw // 2, py + ph + 10), xl, fill="black", font=fax, anchor="ma")

for k in panels:
    draw_panel(k)

out = os.path.join(OUT, "FigS_gating_viability.png")
canvas.save(out, dpi=(300, 300))
print("saved", out, canvas.size)
