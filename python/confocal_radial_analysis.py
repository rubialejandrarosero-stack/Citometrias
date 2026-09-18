#!/usr/bin/env python3
"""Radial distribution of cell populations within spheroids (confocal).

Reproduces the MATLAB DISTRIBUCION.m logic in a single consistent pipeline for the
three channels, from the 2D segmentation masks in
  <BASE>/{time} HORAS/IMAGENES DISTRIBUCION RADIAL/C{channel}_{time}_{sample}.tif
    C1 = PBMC (blue), C2 = MRC-5 (green), C3 = A549 (red), C4 = spheroid (3D stack)
Controls (samples 1 and 4) have no C1 file -> no PBMC, as expected.

For each (time, sample): the spheroid centre & equivalent radius come from C4
(max-projected, Otsu). Each channel mask is labelled; each cell's centroid distance
to the spheroid centre is normalised to the spheroid radius (0 = centre, 100% = edge)
and binned. Output: output/tidy/confocal_radial.csv  (n = 1 spheroid per condition).
"""
import os, numpy as np, csv
from skimage import io, measure, filters

BASE = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/GRAFICAS TESIS/CONFOCAL CELL TRACKER BLUE"
SMAP = {1: ("Basal", "Control"), 2: ("Basal", "Resting"), 3: ("Basal", "Activated"),
        4: ("Inflammatory", "Control"), 5: ("Inflammatory", "Resting"), 6: ("Inflammatory", "Activated")}

def path(c, t, s): return f"{BASE}/{t} HORAS/IMAGENES DISTRIBUCION RADIAL/C{c}_{t}_{s}.tif"
def load(c, t, s):
    p = path(c, t, s); return io.imread(p) if os.path.exists(p) else None
def proj(im): return im.max(0) if (im is not None and im.ndim == 3) else im
def binmask(im):
    if im is None or im.max() == 0: return None
    thr = filters.threshold_otsu(im) if im.max() > im.min() else 0
    return im > thr
def sphere_ref(mask):
    pr = measure.regionprops(measure.label(mask))
    if not pr: return None
    big = max(pr, key=lambda r: r.area); cy, cx = big.centroid
    return cx, cy, np.sqrt(big.area / np.pi)
def centroids_dist(mask, cx, cy, R, minsize=6):
    out = []
    for r in measure.regionprops(measure.label(mask)):
        if r.area < minsize: continue
        cyc, cxc = r.centroid; out.append(np.sqrt((cxc - cx) ** 2 + (cyc - cy) ** 2) / R)
    return np.array(out)

def main():
    bins = np.arange(0, 1.31, 0.1); rows = []
    for t in (24, 48, 96):
        for s in range(1, 7):
            env, cond = SMAP[s]
            ref = sphere_ref(binmask(proj(load(4, t, s))) )
            if ref is None: continue
            for ch, name in [(2, "MRC5"), (3, "A549"), (1, "PBMC")]:
                m = binmask(proj(load(ch, t, s)))
                if m is None: continue            # controls have no C1 -> no PBMC
                h, _ = np.histogram(centroids_dist(m, *ref), bins=bins)
                for i, cnt in enumerate(h):
                    rows.append((t, s, env, cond, name, round((bins[i] + 0.05) * 100), int(cnt)))
    os.makedirs("output/tidy", exist_ok=True)
    with open("output/tidy/confocal_radial.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time", "sample", "environment", "condition", "cell_type", "dist_pct", "count"])
        w.writerows(rows)
    print("Rows:", len(rows))

if __name__ == "__main__":
    main()
