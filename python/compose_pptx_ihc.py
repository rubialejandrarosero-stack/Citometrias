#!/usr/bin/env python3
"""Two IHC compositions from the images in FIGURAS MANUSCRITO.pptx:
  - CD14/CD3  (slide 3): 3 rows (24/48/96 h) x 4 cols [Basal|Inflammatory x resting/activated]
  - PD-L1     (slide 7): 2 rows (Basal/Inflammatory) x 3 cols [Control/resting/activated]
Uniform cells (center-crop to a common aspect), no white balance (preserve DAB). No overlaps.
Outputs: FigS_ihc_CD14_CD3.png (replaces the folder version) / FigS_pdl1.png"""
import io, os
import numpy as np
from scipy.ndimage import gaussian_filter, label, binary_fill_holes, binary_dilation
from pptx import Presentation
from PIL import Image, ImageDraw, ImageFont
from img_norm import white_balance, fit_object


def _largest_cc(m):
    lbl, n = label(m)
    if n == 0:
        return m
    s = np.bincount(lbl.ravel()); s[0] = 0
    return lbl == s.argmax()


def whiten_ihc_bg(im):
    """Segment the stained spheroid (coloured pixels) and replace everything around it — the uneven
    vignette / gradient background, debris and burned-in scale bar — with the UNIFORM pale blue-grey
    slide tone (median of the clean background), so the field reads as one flat colour (not white)."""
    a = np.asarray(im).astype(float)
    sat = a.max(2) - a.min(2)                              # colour saturation (tissue > pale/grey bg)
    val = a.max(2)
    w = im.size[0]
    dens = gaussian_filter((sat > 16).astype(float), sigma=max(4, w // 90))
    body = binary_fill_holes(_largest_cc(dens > 0.22))
    if body.sum() < 500:
        return im
    body = binary_dilation(body, iterations=max(4, w // 55))
    bgmask = (~body) & (val > 200) & (sat < 22)           # clean pale background (no vignette/debris)
    bg = np.median(a[bgmask], axis=0) if bgmask.sum() > 200 else np.array([224.0, 224.0, 230.0])
    out = a.copy(); out[~body] = bg
    return Image.fromarray(out.astype("uint8"))
PPTX = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FIGURAS MANUSCRITO.pptx"
FIG = "output/figuras"
def font(sz):
    return ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", sz)

def slide_pics(idx):
    prs = Presentation(PPTX)
    out = [(round(sh.top / 914400, 1), round(sh.left / 914400, 1), sh.image.blob)
           for sh in prs.slides[idx - 1].shapes if sh.shape_type == 13]
    out.sort(key=lambda x: (x[0], x[1]))
    return [Image.open(io.BytesIO(b)).convert("RGB") for _, _, b in out]

def crop_ar(im, ar):
    """Center-crop to the cell aspect ratio (keeps natural spheroid size)."""
    w, h = im.size
    if w / h > ar:
        nw = int(h * ar); l = (w - nw) // 2; return im.crop((l, 0, l + nw, h))
    nh = int(w / ar); t = (h - nh) // 2
    return im.crop((0, t, w, t + nh))

def build(pics, nrow, ncol, rows, cols, out, title, cw, ch, groups=None, div_after=None, lg=300,
          normalize=False, whiten=False):
    th, gh, ch_h = 66, 58, 56                  # title, group-header, col-header heights
    gap, margin = 16, 44
    grid_w = ncol * cw + (ncol - 1) * gap
    W = margin + lg + grid_w + margin
    y_t = margin
    y_g = y_t + th if groups else y_t + th
    y_c = y_g + (gh if groups else 0)
    y0 = y_c + ch_h
    H = y0 + nrow * ch + (nrow - 1) * gap + margin
    canvas = Image.new("RGB", (W, H), "white"); d = ImageDraw.Draw(canvas)
    def colx(i): return margin + lg + i * (cw + gap)
    d.text((W // 2, y_t), title, fill="black", font=font(56), anchor="ma")
    if groups:
        for lab, c0, c1 in groups:
            xc = (colx(c0) + colx(c1) + cw) // 2
            d.text((xc, y_g), lab, fill="black", font=font(48), anchor="ma")
    for i, c in enumerate(cols):
        d.text((colx(i) + cw // 2, y_c), c, fill="black", font=font(38), anchor="ma")
    for r in range(nrow):
        y = y0 + r * (ch + gap)
        d.text((margin + lg - 16, y + ch // 2), rows[r], fill="black", font=font(46), anchor="rm")
        for c in range(ncol):
            p = pics[r * ncol + c]
            cell = fit_object(p, cw, ch) if normalize else white_balance(crop_ar(p, cw / ch).resize((cw, ch)))
            if whiten:
                cell = whiten_ihc_bg(cell)
            canvas.paste(cell, (colx(c), y))
            # scale bar (~100 µm): black, no number, bottom-right — same style as the H&E montage
            bar = round(0.2145 * cw); bth = max(6, round(cw / 78)); bm = round(cw / 20)
            bx = colx(c) + cw - bm; by = y + ch - bm
            d.rectangle([bx - bar, by - bth, bx, by], fill="black")
    if div_after is not None:
        xd = colx(div_after) + cw + gap // 2
        d.line([(xd, y_c), (xd, y0 + nrow * ch + (nrow - 1) * gap)], fill=(150, 150, 150), width=3)
    canvas.save(os.path.join(FIG, out), dpi=(300, 300))
    print("saved", out, canvas.size)

# --- CD14 / CD3 (slide 3): 3 x 4 ---
build(slide_pics(3), 3, 4,
      rows=["24 h", "48 h", "96 h"],
      cols=["PBMC resting", "PBMC activated", "PBMC resting", "PBMC activated"],
      groups=[("Basal (non-inflammatory)", 0, 1), ("Inflammatory", 2, 3)],
      div_after=None, out="FigS_ihc_CD14_CD3.png",
      title="Immunohistochemistry - CD14 (brown) / CD3 (red)", cw=620, ch=620)

# --- PD-L1 (slide 7): 2 x 3 ---
build(slide_pics(7), 2, 3,
      rows=["Basal", "Inflammatory"],
      cols=["Control (A549/MRC-5)", "PBMC resting", "PBMC activated"],
      out="FigS_pdl1.png", whiten=True,
      title="PD-L1 expression (brown)", cw=820, ch=560, lg=460)

# --- Ki-67 (folder FOTOS IH, 48 h, 6 numbered images): 2 x 3, same layout as PD-L1 ---
IH = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS/FOTOS ESFEROIDES/FOTOS IH"
KI = ["1. Esf no act no inf 48 hrs Ki67.png", "2. Esf +PBMC no act no inf 48 hrs Ki67.png",
      "3. Esf +PBMC act no inf 48 hrs Ki67.png", "4.Esf inf 48 hrs Ki67.png",
      "5. Esf + PBMC no act inf 48 hrs Ki67.png", "6. Esf + PBMC act inf 24 hrs Ki67.png"]  # #6 is 24 h
ki_pics = [Image.open(os.path.join(IH, f)).convert("RGB") for f in KI]
build(ki_pics, 2, 3,
      rows=["Basal", "Inflammatory"],
      cols=["Control (A549/MRC-5)", "PBMC resting", "PBMC activated"],
      out="FigS_ihc_Ki67.png", normalize=True,   # Ki-67: uniform spheroid size on white
      title="Ki-67 expression (brown), 48 h", cw=780, ch=580, lg=460)
