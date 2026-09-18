#!/usr/bin/env python3
"""Gating strategy + population composition over time in one figure:
the gating-strategy panels keep their A-G lettering; the % of CD45+ leukocytes over time is
added below as panel H. Output: output/figuras/FigS_gating_composition.png"""
from PIL import Image, ImageDraw, ImageFont
FIG = "output/figuras"
def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

top = Image.open(f"{FIG}/FigS_gating_strategy.png").convert("RGB")     # panels A-G
bot = Image.open(f"{FIG}/Fig6_unified_por_tiempo.png").convert("RGB")  # panel H

lm, margin, gap = 130, 50, 80
Wc = 3600
tW, tH = Wc, round(Wc * top.height / top.width)
bW, bH = Wc, round(Wc * bot.height / bot.width)
W = lm + Wc + margin
yT = margin; yB = yT + tH + gap
H = yB + bH + margin
canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
canvas.paste(top.resize((tW, tH)), (lm, yT))
canvas.paste(bot.resize((bW, bH)), (lm, yB))
d.text((lm - 96, yB + 6), "H", fill="black", font=font(76), anchor="la")
canvas.save(f"{FIG}/FigS_gating_composition.png", dpi=(300, 300))
print("saved FigS_gating_composition.png", canvas.size)
