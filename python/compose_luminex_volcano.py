#!/usr/bin/env python3
"""Volcano plot of the secretome: Activated vs Resting PBMC (clean contrast — same environment,
differs only in PBMC activation; avoids the exogenous IL-1/TNF confound of Inflammatory vs Basal).
Per cytokine: log2 fold-change (mean Activated / mean Resting) vs -log10(Mann-Whitney p).
Output: output/figuras/Fig_luminex_volcano.png"""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

LAB = {"CCL5 RANTES": "CCL5", "CXCL-9": "CXCL9", "CXCL10/IP10/CRG2": "CXCL10",
       "IL-12/IL-23p40": "IL-12/23", "FGF basic/ FGF2/bFGF": "FGF-2",
       "IL-1beta": "IL-1β", "TNF-alpha": "TNF-α", "IFN-gamma": "IFN-γ", "IFN-alpha": "IFN-α"}
lu = pd.read_csv("output/tidy/luminex_tidy.csv")
lu["grp"] = np.where(lu.pbmc.isin(["Resting", "Activated"]), "PBMC", "Control")

def volcano(gA, gB, get_mask_a, get_mask_b):
    rows = []
    for cyt, g in lu.groupby("cytokine"):
        a = g.loc[get_mask_a(g), "value"].dropna().values
        b = g.loc[get_mask_b(g), "value"].dropna().values
        if len(a) < 3 or len(b) < 3:
            continue
        fc = np.log2((a.mean() + 1) / (b.mean() + 1))
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        rows.append((LAB.get(cyt, cyt), fc, p))
    return pd.DataFrame(rows, columns=["cyt", "log2FC", "p"])

CONTRASTS = {
    "pbmc": ("Secretome: PBMC co-culture vs Control (spheroid only)", "PBMC / Control",
             "Up with PBMC", "Up without PBMC",
             lambda g: g.grp == "PBMC", lambda g: g.grp == "Control"),
    "actrest": ("Secretome: PBMC Activated vs Resting", "Activated / Resting",
                "Up in Activated", "Up in Resting",
                lambda g: g.pbmc == "Activated", lambda g: g.pbmc == "Resting"),
}

def make(key):
    title, xlab, up_lab, dn_lab, ma, mb = CONTRASTS[key]
    df = volcano(None, None, ma, mb)
    df["nlp"] = -np.log10(df.p)
    FC, PS = 1.0, 0.05
    df["sig"] = np.where((df.p < PS) & (df.log2FC >= FC), "up",
                 np.where((df.p < PS) & (df.log2FC <= -FC), "down", "ns"))
    COL = {"up": "#e6550d", "down": "#4292c6", "ns": "#b8b8b8"}
    LEG = {"up": up_lab, "down": dn_lab, "ns": "n.s."}
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    for k, c in COL.items():
        s = df[df.sig == k]
        ax.scatter(s.log2FC, s.nlp, c=c, s=60, edgecolor="black", linewidth=0.4, alpha=0.9, label=LEG[k])
    from adjustText import adjust_text
    txts = [ax.text(r.log2FC, r.nlp, r.cyt, fontsize=9, fontweight="bold")
            for _, r in df[df.sig != "ns"].iterrows()]
    if txts:
        adjust_text(txts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5),
                    expand=(1.3, 1.6))
    ax.axhline(-np.log10(PS), color="grey", ls="--", lw=0.8)
    ax.axvline(FC, color="grey", ls="--", lw=0.8); ax.axvline(-FC, color="grey", ls="--", lw=0.8)
    ax.set_xlabel(f"log2 fold-change ({xlab})", fontsize=13, fontweight="bold")
    ax.set_ylabel("-log10 p (Mann-Whitney)", fontsize=13, fontweight="bold")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(labelsize=11)
    for sp in ["top", "right"]: ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=11, loc="upper left")
    fig.tight_layout()
    out = f"output/figuras/Fig_luminex_volcano_{key}.png"
    fig.savefig(out, dpi=300, facecolor="white"); plt.close(fig)
    print(f"saved {out} | up:{(df.sig=='up').sum()} down:{(df.sig=='down').sum()}")

for key in CONTRASTS:
    make(key)
