#!/usr/bin/env python3
"""Figure 7 - macrophage activation: A = representative CD64/CD11b biaxial density plots,
B = CD64+CD11b+ quantification (Basal), C = CD64+CD11b+ quantification (Inflammatory).
Continuous panel lettering. Output: output/figuras/Fig7_activation.png"""
import os
from PIL import Image, ImageDraw, ImageFont
FIG = "output/figuras"
def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

A = Image.open(f"{FIG}/Fig6_dotplots_representativos.png").convert("RGB")
B = Image.open(f"{FIG}/Fig6act_CD64CD11b_basal.png").convert("RGB")
C = Image.open(f"{FIG}/Fig6act_CD64CD11b_inflammatory.png").convert("RGB")

lm, margin, gap, colgap = 130, 50, 55, 120
Rw = 1450                                   # right-column graph width
Bh = round(Rw * B.height / B.width)
Ch = round(Rw * C.height / C.width)
Ah = Bh + gap + Ch                          # A spans the full height of the stacked B/C column
Aw = round(Ah * A.width / A.height)
W = lm + Aw + colgap + Rw + margin
H = margin + Ah + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
fl = font(76)
# A on the left, full height
canvas.paste(A.resize((Aw, Ah)), (lm, margin))
d.text((lm - 96, margin + 6), "A", fill="black", font=fl, anchor="la")
# right column: B top, C bottom
xR = lm + Aw + colgap
canvas.paste(B.resize((Rw, Bh)), (xR, margin))
canvas.paste(C.resize((Rw, Ch)), (xR, margin + Bh + gap))
d.text((xR - 92, margin + 6), "B", fill="black", font=fl, anchor="la")
d.text((xR - 92, margin + Bh + gap + 6), "C", fill="black", font=fl, anchor="la")
canvas.save(f"{FIG}/Fig7_activation.png", dpi=(300, 300))
print("saved Fig7_activation.png", canvas.size)
