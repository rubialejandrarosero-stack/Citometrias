#!/usr/bin/env python3
"""Secretome + immune-infiltration heatmaps (reference style). Builds two versions:
  full    = all 30 cytokines grouped by family
  key     = only the cytokines significantly induced by PBMC (from the volcano)
Top = cytokine z-score (diverging blue-white-red), bottom = infiltration % of CD45 (yellow-red).
Outputs: output/figuras/Fig_luminex_infiltration.png / Fig_luminex_infiltration_key.png"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

FAM = {
    "Chemokines": ["CCL2", "CCL3", "CCL4", "CCL5 RANTES", "CCL11", "CXCL-9", "CXCL10/IP10/CRG2", "IL-8"],
    "Pro-inflammatory": ["IL-1beta", "IL-6", "TNF-alpha", "IFN-gamma", "IFN-alpha", "IL-12/IL-23p40"],
    "Th2 / regulatory": ["IL-4", "IL-5", "IL-10", "IL-13", "IL-1RA"],
    "T-cell": ["IL-2", "IL-7", "IL-15", "IL-17", "CD25"],
    "Growth factors": ["VEGF", "EGF", "FGF basic/ FGF2/bFGF", "HGF", "G-CSF", "GM-CSF"],
}
# cytokines significantly up with PBMC (volcano): the 'relevant' set
KEY = ["CXCL-9", "CXCL10/IP10/CRG2", "CCL2", "CCL3", "CCL4", "IFN-gamma",
       "IL-10", "IL-1RA", "CD25", "IL-17", "GM-CSF"]
LAB = {"CCL5 RANTES": "CCL5", "CXCL-9": "CXCL9", "CXCL10/IP10/CRG2": "CXCL10",
       "IL-12/IL-23p40": "IL-12/23", "FGF basic/ FGF2/bFGF": "FGF-2",
       "IL-1beta": "IL-1β", "TNF-alpha": "TNF-α", "IFN-gamma": "IFN-γ", "IFN-alpha": "IFN-α"}
COLS = [("Basal", "Control"), ("Basal", "Resting"), ("Basal", "Activated"),
        ("Inflammatory", "Control"), ("Inflammatory", "Resting"), ("Inflammatory", "Activated")]

# ---- matrices ----
lu = pd.read_csv("output/tidy/luminex_tidy.csv")
lu["condition"] = lu["pbmc"].fillna("Control")
M = lu.groupby(["cytokine", "environment", "condition"])["value"].mean().reset_index()
full_mat = M.pivot_table(index="cytokine", columns=["environment", "condition"], values="value").reindex(
    columns=pd.MultiIndex.from_tuples(COLS))

inf = pd.read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv")
inf["environment"] = inf.Condicion.map({"INF": "Inflammatory", "NO_INF": "Basal"})
inf["condition"] = inf.Activacion.map({"ACT": "Activated", "NO_ACT": "Resting"})
POPS = {"pctCD45_CD4": "CD4+ T cells", "pctCD45_CD8": "CD8+ T cells",
        "pctCD45_Mono": "Monocytes (CD14+)", "pctCD45_B": "B cells (CD19+)", "pctCD45_NK": "NK (CD16+)"}
im = inf.groupby(["environment", "condition"])[list(POPS)].mean().rename(columns=POPS)
I = np.full((len(POPS), len(COLS)), np.nan)
for j, (e, c) in enumerate(COLS):
    if (e, c) in im.index:
        I[:, j] = im.loc[(e, c)].values


def build(order, fam, out, fig_h, vmax=1.5):
    mat = full_mat.reindex(order)
    Z = mat.sub(mat.mean(1), axis=0).div(mat.std(1, ddof=1).replace(0, np.nan), axis=0).values
    nC, nI = len(order), len(POPS)
    fig = plt.figure(figsize=(10.4, fig_h))
    gs = GridSpec(2, 1, height_ratios=[nC, nI + 1.2], hspace=0.06, left=0.28, right=0.66,
                  top=0.90, bottom=0.10)
    axc, axi = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    div = axc.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    sc = plt.cm.YlOrRd.copy(); sc.set_bad("#efefef")
    seq = axi.imshow(np.ma.masked_invalid(I), aspect="auto", cmap=sc, vmin=0, vmax=50)
    for ax, rows in [(axc, order), (axi, list(POPS.values()))]:
        ax.set_xticks(range(len(COLS))); ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([LAB.get(r, r) for r in rows], fontsize=9.5, style="italic")
        ax.yaxis.tick_right(); ax.set_xticklabels([]); ax.tick_params(length=0)
        for s in ax.spines.values(): s.set_visible(False)
    for j, (e, c) in enumerate(COLS):
        axc.text(j, -0.9, c, ha="center", va="bottom", fontsize=9, rotation=45)
    axc.text(1, -2.6, "Basal", ha="center", va="bottom", fontsize=12, fontweight="bold")
    axc.text(4, -2.6, "Inflammatory", ha="center", va="bottom", fontsize=12, fontweight="bold")
    axc.axvline(2.5, color="white", lw=3); axi.axvline(2.5, color="white", lw=3)
    if fam:
        row0 = 0
        for f, cyts in fam.items():
            mid = row0 + len(cyts) / 2
            axc.text(-0.10, 1 - mid / nC, f, transform=axc.transAxes, ha="right", va="center",
                     fontsize=10, fontweight="bold", rotation=90)
            if row0 > 0: axc.axhline(row0 - 0.5, color="white", lw=2)
            row0 += len(cyts)
    else:
        axc.text(-0.10, 0.5, "Cytokines", transform=axc.transAxes, ha="right", va="center",
                 fontsize=11, fontweight="bold", rotation=90)
    axi.text(-0.10, 0.5, "Cell\ninfiltration", transform=axi.transAxes, ha="right", va="center",
             fontsize=10, fontweight="bold")
    cax1 = fig.add_axes([0.85, 0.60, 0.020, 0.22]); cb1 = fig.colorbar(div, cax=cax1, ticks=[-vmax, 0, vmax])
    cb1.set_label("z-score", fontsize=10); cb1.ax.tick_params(labelsize=9)
    pi = axi.get_position()                              # anclar la escala al track de infiltración
    cax2 = fig.add_axes([0.85, pi.y0, 0.020, pi.height]); cb2 = fig.colorbar(seq, cax=cax2, ticks=[0, 25, 50])
    cb2.set_label("% of CD45+", fontsize=10); cb2.ax.tick_params(labelsize=9)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("saved", out)


build([c for f in FAM.values() for c in f], FAM, "output/figuras/Fig_luminex_infiltration.png", 12.0)
build(KEY, None, "output/figuras/Fig_luminex_infiltration_key.png", 6.8)
