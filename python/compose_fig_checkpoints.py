#!/usr/bin/env python3
"""Immune-checkpoint figure for the death / immune-blockade experiment, laid out as rows of
[ FlowJo gating layout(s) | its quantification graph ]:
  A  PD-L1+ gate (of CD45- spheroid) | PD-L1+ (% of CD45- spheroid)
  B  PD-1+ gate (of CD3+ T cells)    | CD3+/PD-1+ (% of CD3+ T cells)
Layouts get the same treatment as FigS_gating_viability (crop border+frame, clean axes, headers).
Output: output/figuras/Fig_checkpoints.png"""
import os, numpy as np
from PIL import Image, ImageDraw, ImageFont

LAY = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/LAYOUTS VIABILIDAD MANUSCRITO"
FIG = "output/figuras"
FDIR = "/usr/share/fonts/truetype/dejavu"
def font(sz, bold=True):
    return ImageFont.truetype(os.path.join(FDIR, "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"), sz)

# key: (file, title, parent, x-axis, y-axis)
LAYS = {
    "VIVAS": ("VIVAS", "Live cells",   "all cells",      "FSC-A", "Zombie Red"),
    "CD45":  ("CD45-", "CD45 split",   "live cells",     "FSC-A", "CD45"),
    "PDL1":  ("PD-L1", "PD-L1+ cells", "CD45- spheroid", "FSC-A", "PD-L1"),
    "EPCAM": ("EpCAM", "EpCAM+ tumor", "CD45- spheroid", "FSC-A", "EpCAM"),
    "PD1":   ("PD-1+ ", "PD-1+ T cells", "CD3+ (CD45+)", "CD3",   "PD-1"),
}
# rows: (letter, [layout keys], graph file)
rows = [
    ("A", ["PDL1"], "Fig_pdl1_flow"),
    ("B", ["PD1"],  "Fig_pd1_flow"),
]

def load_plot(name):
    im = Image.open(os.path.join(LAY, name + ".png")).convert("RGBA")
    im = Image.alpha_composite(Image.new("RGBA", im.size, (255, 255, 255, 255)), im).convert("RGB")
    g = np.array(im.convert("L"))
    r = np.where(g.max(axis=1) > 150)[0]; c = np.where(g.max(axis=0) > 150)[0]
    im = im.crop((c.min(), r.min(), c.max() + 1, r.max() + 1)); g = g[r.min():r.max()+1, c.min():c.max()+1]
    H, Wd = g.shape
    vc = [i for i in range(Wd) if (g[:, i] < 80).sum() > 0.6 * H]
    hr = [i for i in range(H) if (g[i, :] < 80).sum() > 0.6 * Wd]
    return im.crop((min(vc), min(hr), max(vc) + 1, max(hr) + 1))

plots = {k: load_plot(v[0]) for k, v in LAYS.items()}

# PD-L1: relabel gate stat (label at left of gate, line at x=173), English + period decimals
_pl = plots["PDL1"]; _dp = ImageDraw.Draw(_pl); _pf = font(28, bold=False)
_dp.rectangle([54, 168, 172, 242], fill="white")
_dp.text((168, 176), "PD-L1+", fill="black", font=_pf, anchor="ra")
_dp.text((168, 208), "9.58", fill="black", font=_pf, anchor="ra")

# PD-1: label sits to the RIGHT of the gate box (right line at x=273); relabel left-aligned there
_p1 = plots["PD1"]; _d1 = ImageDraw.Draw(_p1); _1f = font(28, bold=False)
_d1.rectangle([278, 158, 434, 230], fill="white")
_d1.text((282, 166), "PD-1+", fill="black", font=_1f, anchor="la")
_d1.text((282, 198), "89.7", fill="black", font=_1f, anchor="la")

pw = max(p.width for p in plots.values()); ph = max(p.height for p in plots.values())
ylab_w, header_h, xlab_h = 52, 74, 44
cell_w, cell_h = ylab_w + pw, header_h + ph + xlab_h
fhead, fpar, fax = font(30), font(23, False), font(28)

def layout_cell(key):
    _, title, parent, xl, yl = LAYS[key]
    cell = Image.new("RGB", (cell_w, cell_h), "white"); dd = ImageDraw.Draw(cell)
    dd.text((2, 0), title, fill="black", font=fhead, anchor="la")
    dd.text((2, 38), "from " + parent, fill=(90, 90, 90), font=fpar, anchor="la")
    p = plots[key]; px, py = ylab_w + (pw - p.width) // 2, header_h + (ph - p.height) // 2
    cell.paste(p, (px, py))
    yt = Image.new("RGBA", (ph, ylab_w), (255, 255, 255, 0))
    ImageDraw.Draw(yt).text((ph // 2, ylab_w // 2 + 3), yl, fill="black", font=fax, anchor="mm")
    cell.paste(yt.rotate(90, expand=True), (0, header_h), yt.rotate(90, expand=True))
    dd.text((ylab_w + pw // 2, header_h + ph + 8), xl, fill="black", font=fax, anchor="ma")
    return cell

# Two columns: LEFT = gating layouts (stacked vertical, all same size),
# RIGHT = graphs (stacked vertical, all same size). Rows A/B/C/D read horizontally.
gutter, colgap, rowgap, margin, cgap = 60, 60, 18, 40, 22
GATE_H = 740                                             # gate ~ graph height -> minimal blank between rows
GBOX = (1250, 760)                                       # every graph: same box (fit-to-box)
fl = font(60)

def gate(key):                                           # one gating layout at the uniform size
    c = layout_cell(key)
    return c.resize((round(GATE_H * c.width / c.height), GATE_H))

def gate_block(keys):                                    # stack gate(s) vertically (row A has two)
    gs = [gate(k) for k in keys]
    w = max(g.width for g in gs); h = sum(g.height for g in gs) + cgap * (len(gs) - 1)
    blk = Image.new("RGB", (w, h), "white"); y = 0
    for g in gs:
        blk.paste(g, ((w - g.width) // 2, y)); y += g.height + cgap
    return blk

def graph_fit(gf):                                       # fit graph into the uniform box, centered
    im = Image.open(f"{FIG}/{gf}.png").convert("RGB")
    s = min(GBOX[0] / im.width, GBOX[1] / im.height)
    im = im.resize((round(im.width * s), round(im.height * s)))
    box = Image.new("RGB", GBOX, "white")
    box.paste(im, ((GBOX[0] - im.width) // 2, (GBOX[1] - im.height) // 2))
    return box

built = [(l, gate_block(keys), graph_fit(gf)) for l, keys, gf in rows]
gate_col_w = max(b[1].width for b in built)
xg = margin + gutter                                     # gate column x
xr = xg + gate_col_w + colgap                            # graph column x
W = xr + GBOX[0] + margin
Htot = margin * 2 + sum(max(gb.height, gr.height) for _, gb, gr in built) + rowgap * (len(built) - 1)
canvas = Image.new("RGB", (W, Htot), "white"); d = ImageDraw.Draw(canvas)

y = margin
for letter, gb, gr in built:
    rh = max(gb.height, gr.height)
    gx = xg + (gate_col_w - gb.width) // 2
    gy = y + (rh - gb.height) // 2                       # top of the gating block (dot plots)
    canvas.paste(gb, (gx, gy))
    canvas.paste(gr, (xr, y + (rh - gr.height) // 2))
    d.text((margin + 2, gy - 4), letter, fill="black", font=fl, anchor="la")  # top-left corner
    y += rh + rowgap

out = f"{FIG}/Fig_checkpoints.png"
canvas.save(out, dpi=(300, 300))
print("saved", out, canvas.size, "aspect", round(W / Htot, 2))
