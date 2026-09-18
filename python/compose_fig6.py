#!/usr/bin/env python3
"""Figure 6 - immune populations: A-F = per-population cell-count dot plots
(CD45, CD4, CD8, CD14, NK, B) in a 2x3 grid. The % of CD45+ over time stays a separate
figure (Fig6_unified_por_tiempo). Output: output/figuras/Fig6_populations.png"""
import os
from PIL import Image, ImageDraw, ImageFont
FIG = "output/figuras"
COUNTS = [("A", "Fig6dot_CD45"), ("B", "Fig6dot_CD4"), ("C", "Fig6dot_CD8"),
          ("D", "Fig6dot_Monocytes"), ("E", "Fig6dot_NK"), ("F", "Fig6dot_Bcells")]
def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

colw, gap, lm, margin = 1800, 60, 130, 50
imgs = [(l, Image.open(f"{FIG}/{n}.png").convert("RGB")) for l, n in COUNTS]
ph = round(colw * imgs[0][1].height / imgs[0][1].width)
Wc = 2 * colw + gap
W = lm + Wc + margin
H = margin + 3 * (ph + gap) - gap + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
fl = font(76)
for i, (l, im) in enumerate(imgs):
    r, c = divmod(i, 2)
    x = lm + c * (colw + gap); y = margin + r * (ph + gap)
    canvas.paste(im.resize((colw, ph)), (x, y))
    d.text((x - 96, y + 6), l, fill="black", font=fl, anchor="la")
canvas.save(f"{FIG}/Fig6_populations.png", dpi=(300, 300))
print("saved Fig6_populations.png", canvas.size)
