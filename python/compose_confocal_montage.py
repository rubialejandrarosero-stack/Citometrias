#!/usr/bin/env python3
"""Organise the user's pre-processed confocal montages (from CONFOCAL CLAUDE) into the
manuscript figures, adding unified labels. Each source montage is a 3x4 grid:
rows = Control / Resting / Activated, cols = PBMC(blue) / MRC-5(green) / A549(red) / Merge.

  body : 96 h Non-inflammatory + Inflammatory  -> Fig4_confocal_body.png
  supp : all 6 (24/48/96 h x NI/INF)           -> FigS_confocal_all.png
"""
import os, sys, numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter
from scipy import ndimage


def hq_upscale(im, size):
    """High-quality upscale of a small, heavily-compressed source JPEG: a mild pre-blur removes
    8x8 JPEG block artefacts before LANCZOS resampling, then an unsharp mask restores perceived
    edge sharpness that a plain resize would otherwise leave soft."""
    im = im.filter(ImageFilter.GaussianBlur(radius=0.5))
    im = im.resize(size, resample=Image.LANCZOS)
    return im.filter(ImageFilter.UnsharpMask(radius=2.2, percent=130, threshold=2))

SRC = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/CONFOCAL CLAUDE"
OUT = "output/figuras"
FILES = {
    (24, "NI"): "Montage24 HORAS NO INFLAMATORIO.jpg", (24, "INF"): "Montage24 HORAS INFLAMATORIO.jpg",
    (48, "NI"): "Montage NO INFLAMATORIO 48 HORAS.jpg", (48, "INF"): "Montage INFLAMATORIO 48 HORAS.jpg",
    (96, "NI"): "Montage_96_HORAS_NO_INF.jpg",          (96, "INF"): "Montage_96_HORAS_INF.jpg",
}
COLS = ["PBMC", "MRC-5", "A549", "Merge"]
COLC = [(120, 160, 255), (120, 230, 120), (255, 110, 110), (235, 235, 235)]
ROWS = ["Control", "Resting", "Activated"]
ENVN = {"NI": "Non-inflammatory", "INF": "Inflammatory"}
FDIR = "/usr/share/fonts/truetype/dejavu"
def font(sz, bold=True):
    return ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"), sz)


def clear_green_fog(im, k=58):
    """Raise the black point of the green channel to remove a low-level green background
    haze (saturation glow) so tile backgrounds read as true black. Bright MRC-5 signal is
    preserved; red/blue channels (and their tiles) are untouched."""
    a = np.asarray(im).astype(np.float32)
    a[..., 1] = np.clip((a[..., 1] - k) * (255.0 / (255.0 - k)), 0, 255)
    return Image.fromarray(a.astype("uint8"))


def normalize_spheroids(m, nrows=3, ncols=4, target=0.82):
    """Make every spheroid occupy the same fraction of its cell. For each row the spheroid diameter
    is measured on the Merge cell (last column, always has signal); all 4 channel cells of that row
    are then zoomed by the same factor and re-centred, so rows no longer differ in apparent size."""
    a = np.asarray(m).astype(np.uint8)
    H, W = a.shape[:2]
    cw, ch = W // ncols, H // nrows
    out = a.copy()
    for r in range(nrows):
        mc = a[r * ch:(r + 1) * ch, (ncols - 1) * cw:ncols * cw]
        mask = mc.max(2) > 28
        lbl, n = ndimage.label(mask)
        if n == 0:
            continue
        sizes = np.bincount(lbl.ravel()); sizes[0] = 0
        ys, xs = np.where(lbl == sizes.argmax())
        cx, cy = (xs.min() + xs.max()) / 2, (ys.min() + ys.max()) / 2
        d = max(xs.max() - xs.min(), ys.max() - ys.min()) + 1
        s = target * min(cw, ch) / d
        crop_w, crop_h = cw / s, ch / s
        pad = int(max(crop_w, crop_h)) + 2
        for c in range(ncols):
            cell = a[r * ch:(r + 1) * ch, c * cw:(c + 1) * cw]
            padded = np.zeros((ch + 2 * pad, cw + 2 * pad, 3), np.uint8)
            padded[pad:pad + ch, pad:pad + cw] = cell
            left = int(round(cx + pad - crop_w / 2)); top = int(round(cy + pad - crop_h / 2))
            sub = padded[top:top + int(round(crop_h)), left:left + int(round(crop_w))]
            out[r * ch:(r + 1) * ch, c * cw:(c + 1) * cw] = np.asarray(
                Image.fromarray(sub).resize((cw, ch), resample=Image.LANCZOS))
    return Image.fromarray(out)


def block(time, env, bw):
    """Return a labelled block for one montage. Row labels go in a dedicated left gutter
    (never over the images); the environment tag is rotated on the far left."""
    m = Image.open(os.path.join(SRC, FILES[(time, env)])).convert("RGB")
    if env == "INF":                     # inflammatory montages carry a green background haze
        m = clear_green_fog(m)
    bh = int(bw * m.height / m.width)
    m = normalize_spheroids(hq_upscale(m, (bw, bh)))    # uniform spheroid size across rows
    envw = 46                       # rotated environment tag strip
    roww = int(bw * 0.17)           # row-label gutter (fits "Activated")
    off = envw + roww               # x where the montage starts
    W = off + bw
    canvas = Image.new("RGB", (W, bh), "white"); d = ImageDraw.Draw(canvas)
    canvas.paste(m, (off, 0))
    gw = max(3, bw // 400)              # thin white grid separating every photo (like the H&E montage)
    for k in range(5):                 # 4 columns -> 5 vertical lines (incl. outer border)
        x = min(off + k * bw // 4, off + bw - 1)
        d.line([(x, 0), (x, bh - 1)], fill="white", width=gw)
    for r in range(4):                 # 3 rows -> 4 horizontal lines (incl. outer border)
        y = min(r * bh // 3, bh - 1)
        d.line([(off, y), (off + bw, y)], fill="white", width=gw)
    frow = font(int(bh / 24))
    for r, rl in enumerate(ROWS):
        y = int((r + 0.5) * bh / 3)
        d.text((off - 14, y), rl, fill="black", font=frow, anchor="rm")   # right-aligned in the gutter
    tag = Image.new("RGBA", (bh, envw), (0, 0, 0, 0)); ImageDraw.Draw(tag).text(
        (bh // 2, envw // 2), ENVN[env], fill="black", font=font(int(bh / 20)), anchor="mm")
    canvas.paste(tag.rotate(90, expand=True), (0, 0), tag.rotate(90, expand=True))
    return canvas, off, bw


def col_header(width, off, bw, fsz):
    h = int(fsz * 1.7)
    strip = Image.new("RGB", (width, h), "white"); d = ImageDraw.Draw(strip)
    for c, (name, col) in enumerate(zip(COLS, COLC)):
        x = off + int((c + 0.5) * bw / 4)
        d.text((x, h // 2), name, fill="black", font=font(fsz), anchor="mm")
    return strip


def build(times, out, title=None, bw=1500):
    blocks = []
    off = 0
    for t in times:
        for env in ("NI", "INF"):
            b, off, bwid = block(t, env, bw)
            blocks.append((t, env, b))
    fsz = int(bw / 26)
    hdr = col_header(blocks[0][2].width, off, bw, fsz)
    pad, margin, tgap = 12, 26, int(fsz * 1.5)
    Wc = margin * 2 + blocks[0][2].width
    title_h = int(fsz * 2.0) if title else 0
    time_h = int(fsz * 1.4)
    # layout: title, then per timepoint: time label + col header + NI block + INF block
    per_t = time_h + hdr.height + 2 * blocks[0][2].height + pad
    Hc = margin * 2 + title_h + len(times) * (per_t + tgap)
    canvas = Image.new("RGB", (Wc, Hc), "white"); d = ImageDraw.Draw(canvas)
    y = margin
    if title:
        d.text((Wc // 2, y), title, fill="black", font=font(int(fsz * 1.4)), anchor="ma"); y += title_h
    bi = 0
    for t in times:
        d.text((margin, y), f"{t} h", fill="black", font=font(int(fsz * 1.15)), anchor="la"); y += time_h
        canvas.paste(hdr, (margin, y)); y += hdr.height
        for env in ("NI", "INF"):
            canvas.paste(blocks[bi][2], (margin, y)); y += blocks[bi][2].height; bi += 1
        y += tgap
    canvas.save(os.path.join(OUT, out), dpi=(300, 300))
    print("saved", out, canvas.size)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "body"):
        build([96], "Fig4_confocal_body.png", bw=1600)
    if what in ("all", "supp"):
        for t in (24, 48, 96):
            build([t], f"FigS_confocal_{t}h.png", bw=1600)
