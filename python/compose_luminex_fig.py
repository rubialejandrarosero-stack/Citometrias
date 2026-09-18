#!/usr/bin/env python3
"""Combined secretome figure: A = cytokine heatmap + immune-infiltration track, B = secretome PCA.
Output: output/figuras/Fig_secretome.png"""
from PIL import Image, ImageDraw, ImageFont
FIG = "output/figuras"
def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

A = Image.open(f"{FIG}/Fig_luminex_infiltration.png").convert("RGB")
B = Image.open(f"{FIG}/Fig_luminex_pca.png").convert("RGB")
lm, margin, gap = 120, 50, 90
Bw = 2100; Bh = round(Bw * B.height / B.width); B = B.resize((Bw, Bh))
W = lm + A.width + gap + Bw + margin
H = margin + A.height + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
canvas.paste(A, (lm, margin))
xB = lm + A.width + gap; yB = margin + (A.height - Bh) // 2
canvas.paste(B, (xB, yB))
fl = font(104)
d.text((20, margin + 4), "A", fill="black", font=fl, anchor="la")
d.text((xB - 96, yB + 4), "B", fill="black", font=fl, anchor="la")
canvas.save(f"{FIG}/Fig_secretome.png", dpi=(300, 300))
print("saved Fig_secretome.png", canvas.size)
