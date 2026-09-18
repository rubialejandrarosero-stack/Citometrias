#!/usr/bin/env python3
"""Assemble the representative CD64/CD11b biaxial dot plots (macrophage activation,
gated on CD14+) into a 2x2 grid: rows = environment, columns = PBMC response."""
import os
from PIL import Image, ImageDraw, ImageFont

FIG = "output/figuras"; SRC = os.path.join(FIG, "dotplots")
rows = [("Basal", ["Bas_rest", "Bas_act"]),
        ("Inflammatory", ["Inf_rest", "Inf_act"])]
col_hdr = ["PBMC resting", "PBMC activated"]
pad, bg = 25, "white"
left, top = 300, 150  # margins: row labels (left), column headers + title (top)


def font(sz):
    try:
        return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)
    except OSError:
        return ImageFont.load_default()


tile = Image.open(os.path.join(SRC, f"{rows[0][1][0]}_CD64CD11b.png"))
cw, ch = tile.width, tile.height
W = left + 2 * cw + 3 * pad
H = top + 2 * ch + 3 * pad
canvas = Image.new("RGB", (W, H), bg)
draw = ImageDraw.Draw(canvas)
ftitle, fhdr, frow = font(38), font(34), font(34)

draw.text((left + (2 * cw + pad) // 2 + pad // 2, 34),
          "Macrophage activation — CD64 / CD11b (gated on CD14+)",
          fill="black", font=ftitle, anchor="mm")
for c, htxt in enumerate(col_hdr):
    x = left + pad + c * (cw + pad) + cw // 2
    draw.text((x, top - 34), htxt, fill="black", font=fhdr, anchor="mm")
for r, (rtxt, keys) in enumerate(rows):
    y = top + pad + r * (ch + pad) + ch // 2
    draw.text((left // 2, y), rtxt, fill="black", font=frow, anchor="mm")
    for c, key in enumerate(keys):
        im = Image.open(os.path.join(SRC, f"{key}_CD64CD11b.png")).convert("RGB")
        x = left + pad + c * (cw + pad)
        canvas.paste(im, (x, top + pad + r * (ch + pad)))

out = os.path.join(FIG, "Fig6_dotplots_representativos.png")
canvas.save(out, dpi=(200, 200))
print("Saved", out, canvas.size)
