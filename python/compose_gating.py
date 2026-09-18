#!/usr/bin/env python3
"""Assemble the manuscript gating-strategy figure from the 7 individual FlowJo plot
exports (one representative inflammatory-activated sample). Unifies the style: composites
over white, crops FlowJo's fluorophore axis titles (Comp-*-A) and replaces them with clean
marker names, adds numbered hierarchy headers. HLA-DR plots are intentionally excluded."""
import os
from PIL import Image, ImageDraw, ImageFont

SRC = "output/figuras/gating_dots"          # clean density panels from R/45 (Fig7 style)
OUT = "output/figuras"
# (file, letter, title, parent, x-marker, y-marker)
panels = [
    ("A", "A", "Mononuclear gate",      "all events",     "FSC-A", "SSC-A"),
    ("B", "B", "CD45+ leukocytes",      "mononuclear",    "FSC-A", "CD45"),
    ("C", "C", "T cells (CD3)",         "CD45+",          "FSC-A", "CD3"),
    ("D", "D", "CD4 / CD8 T cells",     "CD3+",           "CD8",   "CD4"),
    ("E", "E", "Monocytes & B cells",   "CD3-",           "CD19",  "CD14"),
    ("F", "F", "NK cells (CD16)",       "CD19- CD14-",    "CD19",  "CD16"),
    ("G", "G", "Macrophage activation", "CD14+",          "CD64",  "CD11b"),
]
FDIR = "/usr/share/fonts/truetype/dejavu"
def font(sz, bold=True):
    return ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"), sz)

def load(name):
    return Image.open(os.path.join(SRC, name + ".png")).convert("RGB")

plots = [load(p[0]) for p in panels]
pw = max(i.width for i in plots); ph = max(i.height for i in plots)

ylab_w, header_h, xlab_h = 54, 78, 48
gx, gy, margin, title_h = 36, 44, 44, 92
cell_w, cell_h = ylab_w + pw, header_h + ph + xlab_h
cols, rows = 4, 2
W = margin * 2 + cols * cell_w + (cols - 1) * gx
H = margin + title_h + rows * cell_h + (rows - 1) * gy + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)

ftitle, fhead, fpar, flet, fax = font(44), font(31), font(24, False), font(36), font(29)
d.text((margin, margin + 20), "Gating strategy", fill="black", font=ftitle, anchor="lm")
d.text((margin, margin + 60), "representative inflammatory, non-activated sample",
       fill=(90, 90, 90), font=fpar, anchor="lm")

top0 = margin + title_h
for k, (name, letter, title, parent, xl, yl) in enumerate(panels):
    r, c = divmod(k, cols)
    cx = margin + c * (cell_w + gx); cy = top0 + r * (cell_h + gy)
    d.text((cx, cy), letter, fill="black", font=flet, anchor="la")
    d.text((cx + 46, cy + 1), title, fill="black", font=fhead, anchor="la")
    d.text((cx + 46, cy + 40), "from " + parent, fill=(90, 90, 90), font=fpar, anchor="la")
    px, py = cx + ylab_w, cy + header_h
    canvas.paste(plots[k], (px, py))
    yt = Image.new("RGBA", (ph, ylab_w), (255, 255, 255, 0))
    ImageDraw.Draw(yt).text((ph // 2, ylab_w // 2 + 4), yl, fill="black", font=fax, anchor="mm")
    yt = yt.rotate(90, expand=True)
    canvas.paste(yt, (cx, py), yt)
    d.text((px + pw // 2, py + ph + 10), xl, fill="black", font=fax, anchor="ma")

out = os.path.join(OUT, "FigS_gating_strategy.png")
canvas.save(out, dpi=(300, 300))
print("saved", out, canvas.size)
