#!/usr/bin/env python3
"""Prepare a figure for Frontiers submission: resize to a target column width at
300 dpi and save as TIFF (LZW). Reusable per figure.

Usage: python make_submission.py <src_basename> <width_mm> [out_basename]
  <src_basename>  PNG in output/figuras/ (without extension)
  <width_mm>      85 (single col), 120 (1.5 col) or 180 (double col)
"""
import sys, os
from PIL import Image

FIG = "output/figuras"; OUT = os.path.join(FIG, "submission")
os.makedirs(OUT, exist_ok=True)
DPI = 300

src = sys.argv[1]
width_mm = float(sys.argv[2])
out = sys.argv[3] if len(sys.argv) > 3 else src

im = Image.open(os.path.join(FIG, f"{src}.png")).convert("RGB")
target_w = round(width_mm / 25.4 * DPI)
target_h = round(im.height * target_w / im.width)
im = im.resize((target_w, target_h), Image.LANCZOS)
dst = os.path.join(OUT, f"{out}.tiff")
im.save(dst, format="TIFF", dpi=(DPI, DPI), compression="tiff_lzw")
mb = os.path.getsize(dst) / 1e6
print(f"Saved {dst}  {target_w}x{target_h}px  = {width_mm}mm x {target_h/DPI*25.4:.0f}mm @ {DPI}dpi  ({mb:.1f} MB)")
