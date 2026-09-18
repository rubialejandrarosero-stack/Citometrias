#!/usr/bin/env python3
"""Figure 1 (morphology, redesigned): left = two brightfield strips (Control vs Activated,
inflammatory) at 24/48/96 h; right = area / diameter / circularity as points joined over time.
Output: output/figuras/Fig1_morphology_full.png"""
import os
from PIL import Image, ImageDraw, ImageFont, ImageFilter
BF = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES/FOTOS ESFEROIDES CONTROL"
FIG = "output/figuras"
# (time label, Control file, Activated file) — inflammatory
ROWS = [
    ("24 h", "24 INF CTL_4.jpg", "24-INF-ACT_6.jpg"),
    ("48 h", "48 INF CTL.jpg",   "48 INF ACT.jpg"),
    ("96 h", "96 INF CTL.jpg",   "96 INF ACT.jpg"),
]
def font(sz, bold=True):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans%s.ttf" % ("-Bold" if bold else ""), sz)

def fit(path, box):
    im = Image.open(path).convert("RGB"); w, h = im.size; s = min(w, h)
    l, t = (w - s) // 2, (h - s) // 2
    return im.crop((l, t, l + s, t + s)).resize((box, box))

cell, gap, rowgap = 700, 16, 16
tgut = 150                                   # left gutter for time labels
envh, condh = 66, 58
margin = 40
grid_w = 2 * cell + gap
left_w = tgut + grid_w
grid_h = 3 * cell + 2 * rowgap
left_h = envh + condh + grid_h               # from top of "Inflammatory" header to bottom row

lines = Image.open(f"{FIG}/preview_morf_lineas.png").convert("RGB")
lines = lines.resize((round(left_h * lines.width / lines.height), left_h))  # match photo-block height

colgap = 70
W = margin + 60 + left_w + colgap + lines.width + margin       # +60 for panel-letter gutter
H = margin + left_h + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)

x0 = margin + 60                              # left block start (after letter gutter)
y_env = margin; y_cond = y_env + envh; y_grid = y_cond + condh
def colx(i): return x0 + tgut + i * (cell + gap)

fenv, fcond, ftime, flet = font(48), font(44), font(40), font(64)
d.text((x0 + tgut + grid_w // 2, y_env), "Inflammatory", fill="black", font=fenv, anchor="ma")
d.text((colx(0) + cell // 2, y_cond + 2), "Control", fill="black", font=fcond, anchor="ma")
d.text((colx(1) + cell // 2, y_cond + 2), "Activated", fill="black", font=fcond, anchor="ma")
fscale = font(24, bold=False)
for r, (tlab, cf, af) in enumerate(ROWS):
    y = y_grid + r * (cell + rowgap)
    cimg = fit(os.path.join(BF, cf), cell); aimg = fit(os.path.join(BF, af), cell)
    for img in (cimg, aimg):                                 # remove burned-in bar: copy bg texture from above, feathered
        bw, bh = 285, 62
        x1, y1 = cell - bw, cell - bh
        src = img.crop((x1, y1 - bh - 8, cell, y1 - 8))
        mask = Image.new("L", (bw, bh), 0)
        ImageDraw.Draw(mask).rectangle([12, 12, bw, bh], fill=255)
        mask = mask.filter(ImageFilter.GaussianBlur(7))
        reg = img.crop((x1, y1, cell, cell)); reg.paste(src, (0, 0), mask)
        img.paste(reg, (x1, y1))
    canvas.paste(cimg, (colx(0), y)); canvas.paste(aimg, (colx(1), y))
    d.text((x0 + tgut - 18, y + cell // 2), tlab, fill="black", font=ftime, anchor="rm")
# single scale bar: exactly the Control-24h calibrated bar (83 px = 200 µm), placed on Activated 96 h
BAR_PX = 83
cbx = colx(1) + cell - 30
cby = y_grid + 2 * (cell + rowgap) + cell - 40
d.rectangle([cbx - BAR_PX, cby - 5, cbx, cby + 5], fill="black")
bar_mid = (cbx - BAR_PX + cbx) / 2                      # exact bar centre
d.text((bar_mid, cby - 11), "200 µm", fill="black", font=fscale, anchor="mb")  # black, no stroke, centred on bar

# right: line-graph panel, vertically centered on the photo grid
xl = x0 + left_w + colgap
yl = y_grid + (grid_h - lines.height) // 2
canvas.paste(lines, (xl, yl))

d.text((10, y_grid - 4), "A", fill="black", font=flet, anchor="lb")
for lab, frac in zip("BCD", (0.21, 0.50, 0.79)):     # one letter per metric row
    d.text((xl - 50, round(yl + frac * lines.height)), lab, fill="black", font=flet, anchor="lm")

canvas.save(f"{FIG}/Fig1_morphology_full.png", dpi=(300, 300))
print("saved Fig1_morphology_full.png", canvas.size, "aspect", round(W / H, 2))
