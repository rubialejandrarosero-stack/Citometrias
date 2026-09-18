#!/usr/bin/env python3
"""PCA of the Luminex secretome, split by environment (Basal / Inflammatory), same style as
Fig_luminex_pca. Each PCA is fit on that environment's samples only (30 cytokines, log10 + z-scored).
Colour = condition, marker = time. Outputs: Fig_luminex_pca_basal.png / _inflammatory.png"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import matplotlib.transforms as tf
from sklearn.decomposition import PCA

plt.rcParams.update({"font.family": "DejaVu Sans"})

lu = pd.read_csv("output/tidy/luminex_tidy.csv")
lu["condition"] = lu["pbmc"].fillna("Control")
lu["rep"] = lu.groupby(["cytokine", "environment", "condition", "time"]).cumcount()
lu["sample"] = (lu.environment + "|" + lu.condition + "|" + lu.time.astype(str) + "|r" + lu.rep.astype(str))

COND_COL = {"Control": "#9e9e9e", "Resting": "#4292c6", "Activated": "#e6550d"}
COND_LAB = {"Control": "Without PBMC", "Resting": "Resting PBMC", "Activated": "Activated PBMC"}
TIME_MRK = {24: "o", 48: "^", 96: "s"}
TIME_LAB = {24: "24 h", 48: "48 h", 96: "96 h"}
ENV_TITLE = {"Basal": "Basal", "Inflammatory": "TNF α, IL-α, IL-1 β"}


def conf_ellipse(ax, x, y, col, n_std=2.0):
    if len(x) < 3: return
    cov = np.cov(x, y); pear = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    rx, ry = np.sqrt(1 + pear), np.sqrt(1 - pear)
    e = Ellipse((0, 0), 2 * rx, 2 * ry, facecolor=col, edgecolor=col, alpha=0.12, lw=1.6, zorder=0)
    sx, sy = np.sqrt(cov[0, 0]) * n_std, np.sqrt(cov[1, 1]) * n_std
    e.set_transform(tf.Affine2D().rotate_deg(45).scale(sx, sy).translate(x.mean(), y.mean()) + ax.transData)
    ax.add_patch(e)


def build(env, out):
    d = lu[lu.environment == env]
    X = d.pivot_table(index="sample", columns="cytokine", values="value")
    meta = d.drop_duplicates("sample").set_index("sample").loc[X.index, ["condition", "time"]]
    X = X.dropna(axis=0); meta = meta.loc[X.index]

    Xz = np.log10(X.clip(lower=0.01))
    Xz = (Xz - Xz.mean()) / Xz.std(ddof=0)
    pc0 = PCA(n_components=2).fit_transform(Xz)
    keep = pc0[:, 0] > -10                      # drop the same type of extreme outlier as the combined PCA
    n_drop = (~keep).sum()
    Xz, meta = Xz[keep], meta[keep]
    pca = PCA(n_components=2).fit(Xz)
    pc = pca.transform(Xz); ev = pca.explained_variance_ratio_ * 100
    print(f"{env}: outliers dropped {n_drop}, n={X.shape[0] - n_drop}")

    fig, ax = plt.subplots(figsize=(13.5, 11.8))
    for cond, col in COND_COL.items():
        m = (meta.condition == cond).values
        conf_ellipse(ax, pc[m, 0], pc[m, 1], col)
    for cond, col in COND_COL.items():
        for t, mk in TIME_MRK.items():
            m = (meta.condition == cond) & (meta.time == t)
            if m.any():
                ax.scatter(pc[m.values, 0], pc[m.values, 1], c=col, marker=mk, s=70,
                           edgecolor="black", linewidth=0.4, alpha=0.85)
    ax.axhline(0, color="grey", lw=0.6, ls=":"); ax.axvline(0, color="grey", lw=0.6, ls=":")
    ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)", fontsize=40, fontweight="bold")
    ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)", fontsize=40, fontweight="bold")
    ax.set_title(f"Secretome PCA — {ENV_TITLE[env]}", fontsize=40, fontweight="bold")
    ax.tick_params(labelsize=36)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)

    leg1 = [Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=12, label=COND_LAB[k])
            for k, c in COND_COL.items()]
    leg2 = [Line2D([0], [0], marker=m, color="w", markerfacecolor="grey", markeredgecolor="black",
                   markersize=12, label=TIME_LAB[k]) for k, m in TIME_MRK.items()]
    l1 = ax.legend(handles=leg1, title="PBMC", loc="upper left", bbox_to_anchor=(1.02, 1.0),
                   frameon=False, fontsize=36, title_fontsize=36)
    ax.add_artist(l1)
    l2 = ax.legend(handles=leg2, title="Time", loc="upper left", bbox_to_anchor=(1.02, 0.55),
                   frameon=False, fontsize=36, title_fontsize=36)
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", bbox_extra_artists=[l1, l2])
    print("saved", out)


build("Basal", "output/figuras/Fig_luminex_pca_basal.png")
build("Inflammatory", "output/figuras/Fig_luminex_pca_inflammatory.png")
