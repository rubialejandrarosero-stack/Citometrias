"""Shared image helpers for the histology/IHC figure montages:
- white_balance: map each channel's bright point to near-white (unify backgrounds).
- fit_object: detect the spheroid (largest tissue blob), crop just that patch and centre it on a
  white canvas sized so the spheroid has the SAME apparent size in every cell. This also removes
  edge vignettes, scanner borders and distant debris. Works for faint spheroids (adaptive threshold)
  and for images with artefact backgrounds (largest connected component)."""
import numpy as np
from scipy.ndimage import gaussian_filter, label
from PIL import Image


def white_balance(im, target=245.0):
    a = np.asarray(im).astype(float)
    ref = np.clip(np.percentile(a.reshape(-1, 3), 90, axis=0), 60, None)
    scale = np.clip(target / ref, 0.75, 1.5)
    return Image.fromarray(np.clip(a * scale, 0, 255).astype("uint8"))


def fit_object(im, cw, ch, fill=0.82, margin=0.06, wb=True):
    if wb:
        im = white_balance(im)
    a = np.asarray(im).astype(float)
    W, H = im.size
    bright = a.mean(2)
    bg = np.percentile(bright, 90)                     # this image's background level
    tissue = bright < bg - 22                          # adaptive: catches faint spheroids too
    tissue = gaussian_filter(tissue.astype(float), 11) > 0.35
    lab, n = label(tissue)                             # keep only the largest blob (the spheroid)
    if n >= 1:
        sizes = np.bincount(lab.ravel()); sizes[0] = 0
        ys, xs = np.where(lab == sizes.argmax())
        y0, y1 = int(np.percentile(ys, 0.5)), int(np.percentile(ys, 99.5))
        x0, x1 = int(np.percentile(xs, 0.5)), int(np.percentile(xs, 99.5))
    else:
        y0, y1, x0, x1 = 0, H, 0, W
    ext = max(y1 - y0, x1 - x0); m = int(ext * margin)
    y0, x0 = max(y0 - m, 0), max(x0 - m, 0)
    y1, x1 = min(y1 + m, H), min(x1 + m, W)
    patch = im.crop((x0, y0, x1, y1))                  # only the spheroid + small margin (no vignette)
    ar = cw / ch
    short = max(patch.width, patch.height) / fill
    ww, wh = (short * ar, short) if ar >= 1 else (short, short / ar)
    ww, wh = int(round(ww)), int(round(wh))
    canvas = Image.new("RGB", (ww, wh), (245, 245, 245))
    canvas.paste(patch, ((ww - patch.width) // 2, (wh - patch.height) // 2))
    return canvas.resize((cw, ch))
