#!/usr/bin/env python3
"""Stack the 3 radial-overlay figures (24/48/96 h) vertically and place a single
legend in the upper-right, top-aligned with the start of the first plot."""
import os
from PIL import Image, ImageChops

FIG = "output/figuras"
panels = ["Fig4_radial_bare_24h", "Fig4_radial_bare_48h", "Fig4_radial_bare_96h"]
pad, bg = 25, "white"

imgs = [Image.open(os.path.join(FIG, f"{n}.png")).convert("RGB") for n in panels]
cw = max(i.width for i in imgs)

# crop the legend to its bounding box
leg = Image.open(os.path.join(FIG, "Fig4_radial_legend.png")).convert("RGB")
diff = ImageChops.difference(leg, Image.new("RGB", leg.size, "white"))
bbox = diff.getbbox()
if bbox:
    leg = leg.crop(bbox)

# right margin holds the legend, fully outside the plots (no overlap)
right_margin = leg.width + 2 * pad
W = cw + right_margin + 2 * pad
H = sum(i.height for i in imgs) + (len(imgs) + 1) * pad
canvas = Image.new("RGB", (W, H), bg)
y = pad
for im in imgs:
    canvas.paste(im, (pad, y))
    y += im.height + pad

# legend in the right margin, top edge ~ where the first panel starts
panel_top = 220          # px from top of the 24 h image to the first plot panel
lx = pad + cw + pad
ly = pad + panel_top
canvas.paste(leg, (lx, ly))

out = os.path.join(FIG, "Fig4_radial_overlay_all.png")
canvas.save(out, dpi=(300, 300))
print("Saved", out, canvas.size)
