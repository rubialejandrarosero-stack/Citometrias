#!/usr/bin/env python3
"""PCA of the Luminex secretome: each sample (environment x condition x time x replicate) as a
point, from the 30 cytokines (log10 + z-scored). Colour = condition, marker = environment.
Output: output/figuras/Fig_luminex_pca.png"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

plt.rcParams.update({"font.family": "DejaVu Sans"})

lu = pd.read_csv("output/tidy/luminex_tidy.csv")
lu["condition"] = lu["pbmc"].fillna("Control")
lu["rep"] = lu.groupby(["cytokine", "environment", "condition", "time"]).cumcount()
lu["sample"] = (lu.environment + "|" + lu.condition + "|" + lu.time.astype(str) + "|r" + lu.rep.astype(str))
X = lu.pivot_table(index="sample", columns="cytokine", values="value")
meta = lu.drop_duplicates("sample").set_index("sample").loc[X.index, ["environment", "condition", "time"]]
X = X.dropna(axis=0)                                    # keep complete samples
meta = meta.loc[X.index]

Xz = np.log10(X.clip(lower=0.01))
Xz = (Xz - Xz.mean()) / Xz.std(ddof=0)
# drop the single extreme outlier, then recompute PCA on the rest (cleaner projection + ellipses)
pc0 = PCA(n_components=2).fit_transform(Xz)
keep = pc0[:, 0] > -10
n_drop = (~keep).sum()
Xz, meta = Xz[keep], meta[keep]
pca = PCA(n_components=2).fit(Xz)
pc = pca.transform(Xz); ev = pca.explained_variance_ratio_ * 100
print(f"outliers dropped: {n_drop}")

COND_COL = {"Control": "#9e9e9e", "Resting": "#4292c6", "Activated": "#e6550d"}
COND_LAB = {"Control": "Without PBMC", "Resting": "Resting PBMC", "Activated": "Activated PBMC"}
ENV_MRK = {"Basal": "o", "Inflammatory": "^"}
ENV_LAB = {"Basal": "Basal", "Inflammatory": "TNF α, IL-α, IL-1 β"}

def conf_ellipse(ax, x, y, col, n_std=2.0):
    from matplotlib.patches import Ellipse
    import matplotlib.transforms as tf
    if len(x) < 3: return
    cov = np.cov(x, y); pear = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    rx, ry = np.sqrt(1 + pear), np.sqrt(1 - pear)
    e = Ellipse((0, 0), 2 * rx, 2 * ry, facecolor=col, edgecolor=col, alpha=0.12, lw=1.6, zorder=0)
    sx, sy = np.sqrt(cov[0, 0]) * n_std, np.sqrt(cov[1, 1]) * n_std
    e.set_transform(tf.Affine2D().rotate_deg(45).scale(sx, sy).translate(x.mean(), y.mean()) + ax.transData)
    ax.add_patch(e)

fig, ax = plt.subplots(figsize=(7.2, 6.4))
for cond, col in COND_COL.items():
    m = (meta.condition == cond).values
    conf_ellipse(ax, pc[m, 0], pc[m, 1], col)
for cond, col in COND_COL.items():
    for env, mk in ENV_MRK.items():
        m = (meta.condition == cond) & (meta.environment == env)
        if m.any():
            ax.scatter(pc[m.values, 0], pc[m.values, 1], c=col, marker=mk, s=70,
                       edgecolor="black", linewidth=0.4, alpha=0.85)
ax.axhline(0, color="grey", lw=0.6, ls=":"); ax.axvline(0, color="grey", lw=0.6, ls=":")
ax.set_xlabel(f"PC1 ({ev[0]:.1f}%)", fontsize=15, fontweight="bold")
ax.set_ylabel(f"PC2 ({ev[1]:.1f}%)", fontsize=15, fontweight="bold")
ax.set_title("Secretome PCA", fontsize=17, fontweight="bold")
ax.tick_params(labelsize=13)
for s in ["top", "right"]: ax.spines[s].set_visible(False)

# legends: condition (colour) + environment (marker)
from matplotlib.lines import Line2D
leg1 = [Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=12, label=COND_LAB[k]) for k, c in COND_COL.items()]
leg2 = [Line2D([0], [0], marker=m, color="w", markerfacecolor="grey", markeredgecolor="black", markersize=12, label=ENV_LAB[k]) for k, m in ENV_MRK.items()]
l1 = ax.legend(handles=leg1, title="PBMC", loc="upper left", bbox_to_anchor=(1.02, 1.0),
               frameon=False, fontsize=13, title_fontsize=13)
ax.add_artist(l1)
ax.legend(handles=leg2, title="Environment", loc="upper left", bbox_to_anchor=(1.02, 0.55),
          frameon=False, fontsize=13, title_fontsize=13)
fig.tight_layout()
fig.savefig("output/figuras/Fig_luminex_pca.png", dpi=300, facecolor="white")
print("saved Fig_luminex_pca.png | PC1 %.1f%% PC2 %.1f%% | n=%d" % (ev[0], ev[1], X.shape[0]))
