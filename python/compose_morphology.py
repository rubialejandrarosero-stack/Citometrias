#!/usr/bin/env python3
"""Stack the 3 morphology panels (area, diameter, circularity) vertically with A/B/C."""
import os
from PIL import Image, ImageDraw, ImageFont

FIG = "output/figuras"
panels = [("Fig1c_area", "A"), ("Fig1c_diameter", "B"), ("Fig1c_circularity", "C")]
pad, bg = 30, "white"
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 110)
except OSError:
    font = ImageFont.load_default()

imgs = [Image.open(os.path.join(FIG, f"{n}.png")).convert("RGB") for n, _ in panels]
cw = max(im.width for im in imgs)
ch = max(im.height for im in imgs)
W = cw + 2 * pad
H = len(imgs) * ch + (len(imgs) + 1) * pad
canvas = Image.new("RGB", (W, H), bg)
draw = ImageDraw.Draw(canvas)
for i, (im, (_, tag)) in enumerate(zip(imgs, panels)):
    y = pad + i * (ch + pad)
    canvas.paste(im, (pad, y))
    draw.text((pad + 12, y + 8), tag, fill="black", font=font)

out = os.path.join(FIG, "Fig1_morphology.png")
canvas.save(out, dpi=(300, 300))
print("Saved", out, canvas.size)
