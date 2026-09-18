#!/usr/bin/env python3
"""Full secretome figure: A = cytokine + infiltration heatmap (left, tall),
B = PCA, C = volcano (PBMC vs Control), D = key-cytokine dot plots (right column).
Output: output/figuras/Fig_secretome_full.png"""
from PIL import Image, ImageDraw, ImageFont
FIG = "output/figuras"
def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

A = Image.open(f"{FIG}/Fig_luminex_infiltration.png").convert("RGB")
B = Image.open(f"{FIG}/Fig_luminex_pca.png").convert("RGB")
C = Image.open(f"{FIG}/Fig_luminex_volcano_pbmc.png").convert("RGB")
D = Image.open(f"{FIG}/Fig_cytokines_key.png").convert("RGB")

lm, margin, g = 120, 50, 70
bw = 1500                                   # B and C width
B = B.resize((bw, round(bw * B.height / B.width)))
C = C.resize((bw, round(bw * C.height / C.width)))
top_h = max(B.height, C.height)
Dw = 2 * bw + g
D = D.resize((Dw, round(Dw * D.height / D.width)))
right_h = top_h + g + D.height
Htot = margin * 2 + max(A.height, right_h)
xR = lm + A.width + g
W = xR + Dw + margin
canvas = Image.new("RGB", (W, Htot), "white"); d = ImageDraw.Draw(canvas)

yA = margin + (max(A.height, right_h) - A.height) // 2
canvas.paste(A, (lm, yA))
canvas.paste(B, (xR, margin))
canvas.paste(C, (xR + bw + g, margin))
yD = margin + top_h + g
canvas.paste(D, (xR, yD))

fl = font(112)
d.text((18, yA + 2), "A", fill="black", font=fl, anchor="la")
d.text((xR - 66, margin + 2), "B", fill="black", font=fl, anchor="la")
d.text((xR + bw + g - 66, margin + 2), "C", fill="black", font=fl, anchor="la")
d.text((xR - 66, yD + 2), "D", fill="black", font=fl, anchor="la")
canvas.save(f"{FIG}/Fig_secretome_full.png", dpi=(300, 300))
print("saved Fig_secretome_full.png", canvas.size)
