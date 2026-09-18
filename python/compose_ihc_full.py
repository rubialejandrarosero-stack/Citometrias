#!/usr/bin/env python3
"""Assemble the complete IHC manuscript figures: panel A = image composition, panel B =
quantification chart, stacked vertically with A/B panel letters. One figure per marker set.
Outputs: Fig_IHC_CD14CD3.png / Fig_IHC_Ki67.png / Fig_IHC_PDL1.png"""
import os
from PIL import Image, ImageDraw, ImageFont
FIG = "output/figuras"
CFG = [  # (image panel, chart panel, output, chart width fraction of content width)
    ("FigS_ihc_Ki67",     "Fig_quant_Ki67",    "Fig_IHC_Ki67.png",    0.52),
    ("FigS_pdl1",         "Fig_quant_PDL1",    "Fig_IHC_PDL1.png",    0.52),
]
def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

# ---- CD14/CD3: images left, Basal chart top-right, Inflammatory chart bottom-right ----
A = Image.open(os.path.join(FIG, "FigS_ihc_CD14_CD3.png")).convert("RGB")
Bg = Image.open(os.path.join(FIG, "Fig_quant_CD14CD3_basal.png")).convert("RGB")
Cg = Image.open(os.path.join(FIG, "Fig_quant_CD14CD3_inflammatory.png")).convert("RGB")
lm, margin, gap, colgap = 120, 50, 60, 110
Rw = 1600
Bh = round(Rw * Bg.height / Bg.width); Ch = round(Rw * Cg.height / Cg.width)
Ah = Bh + gap + Ch
Aw = round(Ah * A.width / A.height)
W = lm + Aw + colgap + Rw + margin
H = margin * 2 + Ah
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
fl = font(96)
canvas.paste(A.resize((Aw, Ah)), (lm, margin))
d.text((22, margin - 6), "A", fill="black", font=fl, anchor="la")
xR = lm + Aw + colgap
canvas.paste(Bg.resize((Rw, Bh)), (xR, margin))
canvas.paste(Cg.resize((Rw, Ch)), (xR, margin + Bh + gap))
d.text((xR - 95, margin - 6), "B", fill="black", font=fl, anchor="la")
d.text((xR - 95, margin + Bh + gap - 6), "C", fill="black", font=fl, anchor="la")
canvas.save(os.path.join(FIG, "Fig_IHC_CD14CD3.png"), dpi=(300, 300))
print("saved Fig_IHC_CD14CD3.png", canvas.size)

for img_name, chart_name, out, wfrac in CFG:
    A = Image.open(os.path.join(FIG, img_name + ".png")).convert("RGB")
    B = Image.open(os.path.join(FIG, chart_name + ".png")).convert("RGB")
    W = A.width                                   # content width = image width
    Bw = int(W * wfrac); Bh = int(B.height * Bw / B.width)
    B = B.resize((Bw, Bh))
    lm, margin, gap = 120, 50, 70                 # left letter margin, outer margin, A/B gap
    Wc = lm + W + margin
    yA = margin; yB = yA + A.height + gap
    Hc = yB + Bh + margin
    canvas = Image.new("RGB", (Wc, Hc), "white"); d = ImageDraw.Draw(canvas)
    canvas.paste(A, (lm, yA))
    canvas.paste(B, (lm + (W - Bw) // 2, yB))
    fl = font(96)
    d.text((22, yA - 6), "A", fill="black", font=fl, anchor="la")
    d.text((22, yB - 6), "B", fill="black", font=fl, anchor="la")
    canvas.save(os.path.join(FIG, out), dpi=(300, 300))
    print("saved", out, canvas.size)
