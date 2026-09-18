"""Shared helpers for the confocal spheroid figures (CellTracker set).
Channels: RED = A549 (tumour), GREEN = MRC-5 (fibroblast), BLUE = PBMC (immune).
Presentation TIFFs live in three folders (24 h and 48 h use *_Z.tif; 96 h drops the _Z)."""
import os, numpy as np
from PIL import Image, ImageSequence
from scipy.ndimage import gaussian_filter
Image.MAX_IMAGE_PIXELS = None

FT = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/FOTOS ESFEROIDES TESIS"
DIRS = {24: f"{FT}/IMAGENES CONFOCAL 24hrs/IMAGENES PRESENTACION24",
        48: f"{FT}/IMAGENES CONFOCAL48 y 96 hrs/IMAGENES PRESENTACION",
        96: f"{FT}/IMAGENES CONFOCAL48 y 96 hrs/CONFOCAL 96 HORAS/IMÁGENES PRESENTACIÓN 96 HORAS"}
# sample number -> (environment, condition)  [from python/confocal_radial_analysis.py SMAP]
SMAP = {1: ("Basal", "Control"), 2: ("Basal", "Resting"), 3: ("Basal", "Activated"),
        4: ("Inflammatory", "Control"), 5: ("Inflammatory", "Resting"), 6: ("Inflammatory", "Activated")}
CH_RGB = {"R": (1, 0, 0), "G": (0, 1, 0), "B": (0, 0, 1)}
CH_NAME = {"R": "RED", "G": "GREEN", "B": "BLUE"}
PIXEL_UM = 1.0 / 0.805152   # ~1.242 µm/px, from ImageJ TIFF metadata (unit=micron)


def crop_bbox(merge, pad=1.18, thr=0.05):
    """Square bounding box (y0,y1,x0,x1) around the spheroid, from the merged intensity."""
    inten = gaussian_filter(merge.sum(2), 8)
    ys, xs = np.where(inten > inten.max() * thr)
    if len(ys) == 0:
        h, w = merge.shape[:2]; return 0, h, 0, w
    cy, cx = (ys.min() + ys.max()) // 2, (xs.min() + xs.max()) // 2
    half = int(max(ys.max() - ys.min(), xs.max() - xs.min()) / 2 * pad)
    h, w = merge.shape[:2]
    half = min(half, cy, cx, h - cy, w - cx)
    return cy - half, cy + half, cx - half, cx + half


def chan_path(time, sample, ch):
    suf = "" if time == 96 else "_Z"
    return os.path.join(DIRS[time], f"{sample}_{CH_NAME[ch]}{suf}.tif")


def _load1(p):
    return np.array(next(ImageSequence.Iterator(Image.open(p)))).astype(float)


def load_corr(time, sample, ch, bg_sigma=60, pct=99.6):
    """Load one channel, subtract a diffuse (gaussian) background, normalise to [0,1]."""
    a = _load1(chan_path(time, sample, ch))
    a = np.clip(a - gaussian_filter(a, bg_sigma), 0, None)
    mx = np.percentile(a, pct)
    return np.clip(a / mx, 0, 1) if mx > 0 else a * 0


def colorize(a, rgb):
    return np.stack([a * c for c in rgb], axis=-1)


def channels_and_merge(time, sample, **kw):
    ch = {c: load_corr(time, sample, c, **kw) for c in "RGB"}
    merge = np.clip(sum(colorize(ch[c], CH_RGB[c]) for c in "RGB"), 0, 1)
    return ch, merge


def u8(arr):
    return (np.clip(arr, 0, 1) * 255).astype("uint8")
