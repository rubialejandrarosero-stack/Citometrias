#!/usr/bin/env python3
"""Unified Figure 4: (A) confocal montage 96 h, (B) cell-type composition ('intensities'),
(C) radial distribution 96 h. A occupies the tall left column; B (top) and C (bottom) stack
on the right. Output: output/figuras/Fig4_confocal_infiltration.png"""
import os
from PIL import Image, ImageDraw, ImageFont
FIG = "output/figuras"
A = Image.open(os.path.join(FIG, "Fig4_confocal_body.png")).convert("RGB")
B = Image.open(os.path.join(FIG, "Fig4_confocal_composition.png")).convert("RGB")
C = Image.open(os.path.join(FIG, "Fig4_radial_body_96h.png")).convert("RGB")

R, g, margin, ml = 2500, 80, 70, 120   # right-col width, gap, outer margin, left letter margin
Bh = round(R * B.height / B.width)
Ch = round(R * C.height / C.width)
Ah = Bh + Ch + g
Aw = round(Ah * A.width / A.height)
A = A.resize((Aw, Ah)); B = B.resize((R, Bh)); C = C.resize((R, Ch))

Wc = ml + Aw + g + R + margin
Hc = margin * 2 + Ah
canvas = Image.new("RGB", (Wc, Hc), "white"); d = ImageDraw.Draw(canvas)
canvas.paste(A, (ml, margin))
xr = ml + Aw + g
canvas.paste(B, (xr, margin))
canvas.paste(C, (xr, margin + Bh + g))

def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)
fl = font(72)
d.text((16, margin + 6), "A", fill="black", font=fl, anchor="la")          # in the left margin, off the image
d.text((xr + 6, margin + 2), "B", fill="black", font=fl, anchor="la")
d.text((xr + 6, margin + Bh + g + 2), "C", fill="black", font=fl, anchor="la")

canvas.save(os.path.join(FIG, "Fig4_confocal_infiltration.png"), dpi=(300, 300))
print("saved Fig4_confocal_infiltration.png", canvas.size)
