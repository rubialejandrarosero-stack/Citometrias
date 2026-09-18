#!/usr/bin/env python3
"""Secretome + immune-infiltration heatmaps, split by environment (same style as
Fig_luminex_infiltration), each showing all 3 timepoints x Control/Resting/Activated (9 columns).
Font +1pt vs the original. Top = cytokine z-score (per row, within this environment's 9 columns),
bottom = infiltration % of CD45. Outputs: Fig_luminex_infiltration_basal.png / _inflammatory.png"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist, squareform

FAM = {
    "Chemokines": ["CCL2", "CCL3", "CCL4", "CCL5 RANTES", "CCL11", "CXCL-9", "CXCL10/IP10/CRG2", "IL-8"],
    "Pro-inflammatory": ["IL-1beta", "IL-6", "TNF-alpha", "IFN-gamma", "IFN-alpha", "IL-12/IL-23p40"],
    "Th2 / regulatory": ["IL-4", "IL-5", "IL-10", "IL-13", "IL-1RA"],
    "T-cell": ["IL-2", "IL-7", "IL-15", "IL-17", "CD25"],
    "Growth factors": ["VEGF", "EGF", "FGF basic/ FGF2/bFGF", "HGF", "G-CSF", "GM-CSF"],
}
LAB = {"CCL5 RANTES": "CCL5", "CXCL-9": "CXCL9", "CXCL10/IP10/CRG2": "CXCL10",
       "IL-12/IL-23p40": "IL-12/23", "FGF basic/ FGF2/bFGF": "FGF-2",
       "IL-1beta": "IL-1β", "TNF-alpha": "TNF-α", "IFN-gamma": "IFN-γ", "IFN-alpha": "IFN-α"}
TIMES = [24, 48, 96]
CONDS = ["Control", "Resting", "Activated"]
COLS = [(t, c) for t in TIMES for c in CONDS]                 # 9 columns: time-major, cond-minor
ORDER = [c for f in FAM.values() for c in f]

lu = pd.read_csv("output/tidy/luminex_tidy.csv")
lu["condition"] = lu["pbmc"].fillna("Control")

inf = pd.read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv")
inf["environment"] = inf.Condicion.map({"INF": "Inflammatory", "NO_INF": "Basal"})
inf["condition"] = inf.Activacion.map({"ACT": "Activated", "NO_ACT": "Resting"})
POPS = {"pctCD45_CD4": "CD4+ T cells", "pctCD45_CD8": "CD8+ T cells",
        "pctCD45_Mono": "CD14+ cells", "pctCD45_B": "CD19+ cells", "pctCD45_NK": "CD16+ cells"}


def permanova(D, groups, n_perm=9999, seed=0):
    """PERMANOVA (Anderson 2001) pseudo-F on a square distance matrix D for a single factor
    `groups`, with permutation-based p-value. Returns (R2, F, p)."""
    rng = np.random.default_rng(seed)
    D = np.asarray(D)
    n = D.shape[0]
    D2 = D ** 2
    ss_total = D2.sum() / (2 * n)

    def ss_within(order):
        g = np.asarray(groups)[order]
        ssw = 0.0
        for lev in np.unique(g):
            idx = np.where(g == lev)[0]
            if len(idx) > 1:
                ssw += D2[np.ix_(idx, idx)].sum() / (2 * len(idx))
        return ssw

    idx0 = np.arange(n)
    ssw0 = ss_within(idx0)
    ssa0 = ss_total - ssw0
    a = len(np.unique(groups)) - 1
    dfw = n - a - 1
    F0 = (ssa0 / a) / (ssw0 / dfw)
    R2 = ssa0 / ss_total

    ge = 1
    for _ in range(n_perm):
        perm = rng.permutation(n)
        ssw = ss_within(perm)
        ssa = ss_total - ssw
        F = (ssa / a) / (ssw / dfw)
        if F >= F0:
            ge += 1
    p = ge / (n_perm + 1)
    return R2, F0, p


def run_permanova_env(lu, env):
    """PERMANOVA on Euclidean distance of the z-scored cytokine profile, restricted to Resting +
    Activated (the only conditions with real biological replicates -- 4 donors -- so within-group
    variance is estimable). Control is excluded from the test (kept in the heatmap display only).
    Two single-factor models are reported: condition (Resting vs Activated) and time (24/48/96h)."""
    sub = lu[(lu.environment == env) & (lu.pbmc.isin(["Resting", "Activated"]))]
    wide = sub.pivot_table(index=["condition", "time", "donor"], columns="cytokine", values="value")
    wide = wide.dropna(axis=1, how="any")  # complete cases only, so distances are comparable
    Zs = (wide - wide.mean()) / wide.std(ddof=1)
    D = squareform(pdist(Zs.values, metric="euclidean"))
    meta = wide.index.to_frame(index=False)
    res = {}
    for factor in ("condition", "time"):
        R2, F, p = permanova(D, meta[factor].values)
        res[factor] = (R2, F, p)
    return res, wide.shape[0]


def cluster_within_family(mat_z, fam):
    """Reorder cytokines within each family block by hierarchical clustering (average linkage,
    correlation distance) on their z-score profile across the 9 columns, so families/labels are
    preserved but rows within a block follow similarity of pattern rather than input order."""
    order = []
    for cyts in fam.values():
        sub = mat_z.loc[cyts]
        if len(cyts) <= 2 or sub.isna().any(axis=1).any():
            order.extend(cyts)  # too few rows, or missing values -> clustering undefined, keep as-is
            continue
        d = pdist(sub.values, metric="correlation")
        if not np.isfinite(d).all():
            order.extend(cyts)
            continue
        leaves = leaves_list(linkage(d, method="average"))
        order.extend([cyts[i] for i in leaves])
    return order


def build(env, out, vmax=1.5):
    perm_res, n_perm = run_permanova_env(lu, env)
    Mlu = lu[lu.environment == env].groupby(["cytokine", "time", "condition"])["value"].mean().reset_index()
    mat0 = Mlu.pivot_table(index="cytokine", columns=["time", "condition"], values="value").reindex(
        index=ORDER, columns=pd.MultiIndex.from_tuples(COLS))
    z0 = mat0.sub(mat0.mean(1), axis=0).div(mat0.std(1, ddof=1).replace(0, np.nan), axis=0)
    order = cluster_within_family(z0, FAM)
    mat = mat0.loc[order]
    Z = z0.loc[order].values

    ie = inf[inf.environment == env].groupby(["Tiempo", "condition"])[list(POPS)].mean().rename(columns=POPS)
    I = np.full((len(POPS), len(COLS)), np.nan)
    for j, (t, c) in enumerate(COLS):
        if (t, c) in ie.index:
            I[:, j] = ie.loc[(t, c)].values

    nC, nI = len(ORDER), len(POPS)
    fig = plt.figure(figsize=(32.0, 40.0))
    gs = GridSpec(2, 1, height_ratios=[nC, nI + 1.2], hspace=0.06, left=0.27, right=0.62,
                  top=0.84, bottom=0.06)
    axc, axi = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    div = axc.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    sc = plt.cm.YlOrRd.copy(); sc.set_bad("#efefef")
    seq = axi.imshow(np.ma.masked_invalid(I), aspect="auto", cmap=sc, vmin=0, vmax=50)
    for ax, rows in [(axc, order), (axi, list(POPS.values()))]:
        ax.set_xticks(range(len(COLS))); ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([LAB.get(r, r) for r in rows], fontsize=40, style="italic")
        ax.yaxis.tick_right(); ax.set_xticklabels([]); ax.tick_params(length=0)
        for s in ax.spines.values(): s.set_visible(False)
    cond_lab = {"Control": "Control", "Resting": "Resting PBMC", "Activated": "Activated PBMC"}
    COND_COL = {"Control": "black", "Resting": "#4292c6", "Activated": "#e31a1c"}
    for j, (t, c) in enumerate(COLS):
        axc.add_patch(plt.Rectangle((j - 0.5, -1.5), 1, 1, transform=axc.transData,
                                     facecolor=COND_COL[c], edgecolor="white", linewidth=1.5,
                                     clip_on=False))
    for k, t in enumerate(TIMES):
        axc.text(1 + k * 3, -3.0, f"{t} h", ha="center", va="bottom", fontsize=40, fontweight="bold")
    env_lab = env
    axc.text(4, -5.1, env_lab, ha="center", va="bottom", fontsize=40, fontweight="bold")
    for x in (2.5, 5.5):
        axc.axvline(x, color="white", lw=3); axi.axvline(x, color="white", lw=3)
    FAM_LAB = {"Th2 / regulatory": "Th2-\nregulatory", "Pro-inflammatory": "Pro-\ninflammatory",
               "Growth factors": "Growth\nfactors"}
    row0 = 0
    for f, cyts in FAM.items():
        mid = row0 + len(cyts) / 2
        flab = FAM_LAB.get(f, f)
        axc.text(-0.11, 1 - mid / nC, flab, transform=axc.transAxes, ha="right", va="center",
                 fontsize=32, fontweight="bold", rotation=0, linespacing=1.1)
        if row0 > 0: axc.axhline(row0 - 0.5, color="white", lw=2)
        row0 += len(cyts)
    axi.text(-0.11, 0.5, "Cell\ninfiltration", transform=axi.transAxes, ha="right", va="center",
             fontsize=40, fontweight="bold")
    cax1 = fig.add_axes([0.92, 0.60, 0.020, 0.22]); cb1 = fig.colorbar(div, cax=cax1, ticks=[-vmax, 0, vmax])
    cb1.set_label("z-score", fontsize=40); cb1.ax.tick_params(labelsize=36)
    pi = axi.get_position()
    cb2_h = pi.height * 0.85
    cax2 = fig.add_axes([0.92, pi.y0 + (pi.height - cb2_h) / 2, 0.020, cb2_h])
    cb2 = fig.colorbar(seq, cax=cax2, ticks=[0, 25, 50])
    cb2.set_label("% of CD45+", fontsize=40, labelpad=20); cb2.ax.tick_params(labelsize=36)

    leg_ax = fig.add_axes([0.885, cax1.get_position().y0 - 0.115, 0.11, 0.09])
    leg_ax.axis("off")
    for k, cnd in enumerate(CONDS):
        y = 1 - k / (len(CONDS) - 1)
        leg_ax.add_patch(plt.Rectangle((0.0, y - 0.07), 0.16, 0.14, transform=leg_ax.transAxes,
                                        facecolor=COND_COL[cnd], edgecolor="white", linewidth=1.0,
                                        clip_on=False))
        leg_ax.text(0.24, y, cond_lab[cnd], transform=leg_ax.transAxes, ha="left", va="center",
                    fontsize=30)

    def fmt_p(p):
        return "< 0.0001" if p < 0.0001 else f"= {p:.4f}"
    perm_lines = [f"PERMANOVA (Euclidean, z-score, n={n_perm})", "Resting vs Activated only:"]
    for factor, label in [("condition", "Condition"), ("time", "Time")]:
        R2, F, p = perm_res[factor]
        perm_lines.append(f"{label}: R²={R2:.3f}, F={F:.2f}, p {fmt_p(p)}")
    perm_ax = fig.add_axes([0.885, leg_ax.get_position().y0 - 0.145, 0.11, 0.115])
    perm_ax.axis("off")
    perm_ax.text(0.0, 1.0, "\n".join(perm_lines), transform=perm_ax.transAxes, ha="left", va="top",
                 fontsize=24, linespacing=1.6,
                 bbox=dict(boxstyle="round,pad=0.4", facecolor="#f5f5f5", edgecolor="grey"))

    fig.savefig(out, dpi=220, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("saved", out)


build("Basal", "output/figuras/Fig_luminex_infiltration_basal.png")
build("Inflammatory", "output/figuras/Fig_luminex_infiltration_inflammatory.png")
