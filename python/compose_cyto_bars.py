#!/usr/bin/env python3
"""Assemble the cytometry dot-plot panels (unchanged) into one montage.
CD45 is shown separately (its own figure); this composition has CD4/CD8 on top and
CD14/NK/B cells on the bottom row."""
import os
from PIL import Image, ImageDraw, ImageFont

FIG = "output/figuras"
# rows of (filename, panel-tag); rows may have different numbers of panels
rows = [
    [("Fig6dot_CD4", "A"), ("Fig6dot_CD8", "B")],
    [("Fig6dot_Monocytes", "C"), ("Fig6dot_NK", "D")],
    [("Fig6dot_Bcells", "E")],
]
pad, bg = 40, "white"

imgs = {n: Image.open(os.path.join(FIG, f"{n}.png")).convert("RGB") for row in rows for n, _ in row}
cw = max(im.width for im in imgs.values())
ch = max(im.height for im in imgs.values())
maxcols = max(len(r) for r in rows)
W = maxcols * cw + (maxcols + 1) * pad
H = len(rows) * ch + (len(rows) + 1) * pad

canvas = Image.new("RGB", (W, H), bg)
draw = ImageDraw.Draw(canvas)
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 90)
except OSError:
    font = ImageFont.load_default()

for r, row in enumerate(rows):
    k = len(row)
    x0 = pad  # left-align rows that have fewer panels than the widest row
    y = pad + r * (ch + pad)
    for i, (name, tag) in enumerate(row):
        x = x0 + i * (cw + pad)
        canvas.paste(imgs[name], (x, y))
        draw.text((x + 10, y + 5), tag, fill="black", font=font)

out = os.path.join(FIG, "Fig6_composition_counts.png")
canvas.save(out, dpi=(300, 300))
print("Saved", out, canvas.size)
