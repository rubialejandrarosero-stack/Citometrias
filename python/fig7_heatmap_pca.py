#!/usr/bin/env python3
"""Figura combinada heatmap (30 citoquinas) + PCA, por ambiente (Basal / Inflammatory), con
PERMANOVA de 2 factores (condicion PBMC x tiempo) sobre distancia euclidea (log10 + z-score),
adaptado del diseno de fig7cens.py (Scripts_Figura7) a nuestra fuente de datos verificada
(los dos CSV crudos Luminex xPONENT), en vez de los xlsx intermedios que no tenemos.

Control (Sph/esferoide sin PBMC) entra al PCA y al PERMANOVA con sus 4 replicas tecnicas de
pozo como observaciones individuales (igual que el script original trata cada xlcol), no como
un unico valor promediado -- replica fielmente ese diseno, a diferencia del heatmap
Fig_luminex_infiltration_*, donde Control se promedia y el PERMANOVA excluye Control.
"""
import sys
sys.path.insert(0, "python")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Verdana", "DejaVu Sans"]
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA
from permanova import adonis2
from build_luminex_tidy_from_raw import extract_result_block, RAW_FILES

FAM = [
    ("Chemokines", ["CCL2", "CCL3", "CCL4", "CCL5 RANTES", "CCL11", "CXCL-9", "CXCL10/IP10/CRG2", "IL-8"]),
    ("Pro-inflammatory", ["IL-1beta", "IL-6", "TNF-alpha", "IFN-gamma", "IFN-alpha", "IL-12/IL-23p40"]),
    ("Th2 /\nregulatory", ["IL-4", "IL-5", "IL-10", "IL-13", "IL-1RA"]),
    ("T cell", ["IL-2", "IL-7", "IL-15", "IL-17", "CD25"]),
    ("Growth factors", ["VEGF", "EGF", "FGF basic/ FGF2/bFGF", "HGF", "G-CSF", "GM-CSF"]),
]
ALL = [a for _, g in FAM for a in g]
assert len(ALL) == 30
NM = {"CCL5 RANTES": "CCL5", "CXCL-9": "CXCL9", "CXCL10/IP10/CRG2": "CXCL10", "IL-1beta": "IL-1β",
      "TNF-alpha": "TNF-α", "IFN-gamma": "IFN-γ", "IFN-alpha": "IFN-α",
      "FGF basic/ FGF2/bFGF": "FGF-2", "IL-12/IL-23p40": "IL-12/23"}
ARMS = ["Control", "Resting", "Activated"]
COL = {"Control": "black", "Resting": "#2c7fb8", "Activated": "#d7301f"}
MK = {24: "o", 48: "s", 96: "^"}

# ---- raw long-format data, WITHOUT collapsing Control's 4 technical replicates ----
frames = [extract_result_block(p) for p in RAW_FILES]
X = pd.concat(frames, ignore_index=True)
X["condition"] = X.pbmc.map({"Resting": "Resting", "Activated": "Activated"}).fillna("Control")
# pseudo-well id: for Control, the 4 technical replicates sharing (cytokine, env, time, sample)
# in file order; for PBMC, donor already uniquely identifies the replicate.
X["rep"] = X.groupby(["cytokine", "environment", "time", "sample"]).cumcount() + 1
X["obs_id"] = np.where(X.condition == "Control",
                        X.environment + "_" + X.time.astype(str) + "_C" + X.rep.astype(str),
                        X.environment + "_" + X.time.astype(str) + "_" + X.condition.str[:1] + "_D" + X.donor.astype("Int64").astype(str))

cols = [(env, t, a) for env in ["Basal", "Inflammatory"] for t in [24, 48, 96] for a in ARMS]
M = np.array([[X[(X.cytokine == an) & (X.environment == f) & (X.time == t) & (X.condition == a)].value.mean()
               for f, t, a in cols] for an in ALL])
L10 = np.log10(M + 1)
Z = (L10 - L10.mean(1, keepdims=True)) / L10.std(1, keepdims=True)

fig = plt.figure(figsize=(18, 20.5))
gs = GridSpec(2, 6, figure=fig, height_ratios=[1.62, 1.0], hspace=.30, wspace=.60)
axA = fig.add_subplot(gs[0, :])
im = axA.imshow(Z, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto")
lbl = [NM.get(a, a) for a in ALL]
axA.set_yticks(range(30)); axA.set_yticklabels(lbl, fontsize=15.5, style="italic")
XLAB = {"Control": "Control", "Resting": "Resting PBMC", "Activated": "Activated PBMC"}
axA.set_xticks(range(18))
axA.set_xticklabels([XLAB[a] for _, _, a in cols], fontsize=13, rotation=90)
axA.set_xticks(np.arange(-.5, 18, 1), minor=True); axA.set_yticks(np.arange(-.5, 30, 1), minor=True)
axA.grid(which="minor", color="white", lw=.45); axA.tick_params(which="minor", length=0)
axA.axvline(8.5, c="k", lw=2.8)
for x in [2.5, 5.5, 11.5, 14.5]:
    axA.axvline(x, c="white", lw=2.2)
y0 = 0
for fe, g in FAM:
    y1 = y0 + len(g)
    if y1 < 30:
        axA.axhline(y1 - .5, c="white", lw=2.6)
    axA.annotate(fe, xy=(0, (y0 + y1 - 1) / 2), xycoords=("axes fraction", "data"),
                 xytext=(-122, 0), textcoords="offset points", ha="right", va="center",
                 fontsize=16, fontweight="bold", color="0.3", annotation_clip=False)
    y0 = y1
for k, t in enumerate([24, 48, 96] * 2):
    axA.text(k * 3 + 1, -1.3, f"{t} h", ha="center", fontsize=17, fontweight="bold")
axA.text(4, -2.6, "Basal", ha="center", fontsize=20, fontweight="bold")
axA.text(13, -2.6, "Inflammatory", ha="center", fontsize=20, fontweight="bold")
cb = fig.colorbar(im, ax=axA, fraction=.015, pad=.055); cb.set_label("z-score", fontsize=16)
cb.ax.tick_params(labelsize=14)
axA.text(-.16, .99, "A", transform=axA.transAxes, fontsize=32, fontweight="bold")

RES = {}
for k, (f, ttl) in enumerate([("Basal", "Basal"), ("Inflammatory", "Inflammatory")]):
    s = X[X.environment == f]
    P = s.pivot_table(index="obs_id", columns="cytokine", values="value")[ALL]
    me = s.groupby("obs_id")[["time", "condition"]].first().loc[P.index]
    Y = np.log10(P + 1); Y = (Y - Y.mean()) / Y.std()
    pc = PCA(); S = pc.fit_transform(Y.values); ev = pc.explained_variance_ratio_ * 100
    r = adonis2(Y.values, me.condition.values, me.time.values.astype(str))
    RES[f] = dict(ev=ev, **r, n=len(P))
    ax = fig.add_subplot(gs[1, k * 3:(k + 1) * 3])
    for a in ARMS:
        for t in [24, 48, 96]:
            j = ((me.condition == a) & (me.time == t)).values
            if j.sum():
                ax.scatter(S[j, 0], S[j, 1], c=COL[a], marker=MK[t], s=115, edgecolor="k", lw=.6, zorder=3)
        j = (me.condition == a).values
        ax.scatter(S[j, 0].mean(), S[j, 1].mean(), c=COL[a], marker="X", s=340, edgecolor="k", lw=1.4, zorder=4)
    ax.set_xlabel(f"PC1 ({ev[0]:.1f} %)", fontsize=18); ax.set_ylabel(f"PC2 ({ev[1]:.1f} %)", fontsize=18)
    ax.set_title(ttl, fontweight="bold", fontsize=20.5, pad=10)
    ax.tick_params(labelsize=13)
    ax.axhline(0, lw=.5, c="k"); ax.axvline(0, lw=.5, c="k"); ax.grid(ls=":", alpha=.3)
    a0, a1 = ax.get_ylim(); ax.set_ylim(a0 - (a1 - a0) * .32, a1)
    pf = lambda p: "p < 0.0001" if p < 1e-4 else f"p = {p:.4f}"
    r2 = lambda v: f"R² = {v:.3f}"
    lab = (f"PERMANOVA\nPBMC condition   {r2(r['R2_f1'])}   {pf(r['p_f1'])}"
           f"\ntime                      {r2(r['R2_f2'])}   {pf(r['p_f2'])}")
    ax.text(.03, .03, lab, transform=ax.transAxes, va="bottom", fontsize=14.5,
            bbox=dict(fc="white", ec=".65", alpha=.93, pad=4))
    ax.text(-.135, 1.09, "B" if k == 0 else "C", transform=ax.transAxes, fontsize=32, fontweight="bold")

LBL = {"Control": "Control", "Resting": "Resting PBMC", "Activated": "Activated PBMC"}
H = ([Line2D([], [], marker="o", ls="", color=COL[a], ms=14, label=LBL[a]) for a in ARMS] +
     [Line2D([], [], marker=MK[t], ls="", mfc=".45", mec="k", ms=13, label=f"{t} h") for t in [24, 48, 96]] +
     [Line2D([], [], marker="X", ls="", mfc="w", mec="k", ms=15, label="centroid")])
fig.legend(handles=H, loc="lower center", ncol=7, frameon=False, fontsize=18, bbox_to_anchor=(.5, .006))
plt.subplots_adjust(left=.24, right=.84, top=.955, bottom=.075)

out = "output/figuras/Fig7_heatmap_pca_combined.png"
fig.savefig(out, dpi=220, facecolor="white", bbox_inches="tight")
print("saved", out)
for f in RES:
    r = RES[f]
    print(f"  {f}: n={r['n']} PC1={r['ev'][0]:.1f} PC2={r['ev'][1]:.1f} "
          f"cond R2={r['R2_f1']:.3f} p={r['p_f1']:.4f} | tiempo R2={r['R2_f2']:.3f} p={r['p_f2']:.4f}")
