# R/45_gating_dotplots.R — gating-strategy panels in the Fig7_activation style (pseudocolour density
# from raw FCS), with each population enclosed in a rectangle (bounding box of the cluster) + %.
# Sample: Inflam_Noact_3 (inflammatory, non-activated, 24 h) from the FlowJo .acs archive.
# Panel A uses an FSC threshold only. Output: output/figuras/gating_dots/{A..G}.png
suppressPackageStartupMessages({library(flowCore); library(ggplot2)})

FCS <- "data/raw/acs_24h/Inflam_Noact_3.fcs"
UNS <- "data/raw/Infiltracion/Reference Group/Unstained (Cells).fcs"
OUT <- "output/figuras/gating_dots"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ch <- c(CD45="APC-Cy7-A", CD3="PE-Cy7-A", CD4="BV421-A", CD8="BV510-A", CD14="PE-A",
        CD16="FITC-A", CD19="StarBrightViolet 475-A", CD64="PerCP-Cy5.5-A", CD11b="APC-A")
FSC_MIN <- 5e5

rd  <- function(f) read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
lg  <- estimateLogicle(rd(FCS), channels = unname(ch))
tdf <- function(f) as.data.frame(exprs(transform(rd(f), lg)))
e   <- tdf(FCS)
un  <- tdf(UNS)
th  <- setNames(sapply(ch, function(c) as.numeric(quantile(un[[c]], 0.999, na.rm = TRUE))), names(ch))

flow_ramp <- colorRampPalette(c("#000091", "#0000FF", "#00C0FF", "#00FF60",
                                "#B0FF00", "#FFD000", "#FF3000", "#B00000"))
pc1 <- function(n, d) round(n / max(d, 1) * 100, 1)

base_plot <- function(dat, xcol, ycol) {
  d <- data.frame(x = dat[[xcol]], y = dat[[ycol]])
  d$col <- densCols(d$x, d$y, colramp = flow_ramp, nbin = 256)
  qx <- quantile(d$x, c(0.001, 0.999), na.rm = TRUE); qy <- quantile(d$y, c(0.001, 0.999), na.rm = TRUE)
  list(p = ggplot(d, aes(x, y)) + geom_point(colour = d$col, size = 0.11, shape = 16) +
         coord_cartesian(xlim = qx, ylim = qy, expand = FALSE) +
         labs(x = NULL, y = NULL) + theme_classic(base_size = 15) +
         theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
               plot.margin = margin(3, 3, 3, 3),
               panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6)),
       xl = qx, yl = qy)
}
# bounding box (rectangle) of a sub-population, clamped to the plot window
bbox <- function(dat, xcol, ycol, mask, xl, yl, q = 0.02) {
  s <- dat[mask, , drop = FALSE]
  bx <- quantile(s[[xcol]], c(q, 1 - q), na.rm = TRUE); by <- quantile(s[[ycol]], c(q, 1 - q), na.rm = TRUE)
  c(max(bx[1], xl[1]), min(bx[2], xl[2]), max(by[1], yl[1]), min(by[2], yl[2]))
}
GATE <- function(p, b) p + annotate("rect", xmin = b[1], xmax = b[2], ymin = b[3], ymax = b[4],
                                     fill = NA, colour = "black", linewidth = 0.6)
LAB <- function(p, x, y, txt, hj = 0, vj = 1) p +
  annotate("text", x = x, y = y, label = txt, hjust = hj, vjust = vj, fontface = "bold", size = 4)
# labels at fixed plot corners (avoid clipping): "tl" top-left, "br" bottom-right
LTL <- function(p, xl, yl, txt) LAB(p, xl[1] + 0.03 * diff(xl), yl[2] - 0.03 * diff(yl), txt, hj = 0, vj = 1)
LBR <- function(p, xl, yl, txt) LAB(p, xl[2] - 0.03 * diff(xl), yl[1] + 0.03 * diff(yl), txt, hj = 1, vj = 0)
sv <- function(p, f) ggsave(file.path(OUT, f), p, width = 3.4, height = 3.4, dpi = 200, bg = "white")

# ---- gating hierarchy (threshold-based) ----
mono <- e[e[["FSC-A"]] > FSC_MIN, , drop = FALSE]
cd45 <- mono[mono[[ch["CD45"]]] > th["CD45"], , drop = FALSE]
cd3p <- cd45[cd45[[ch["CD3"]]] > th["CD3"], , drop = FALSE]
cd3n <- cd45[cd45[[ch["CD3"]]] <= th["CD3"], , drop = FALSE]
cd14 <- cd3n[cd3n[[ch["CD14"]]] > th["CD14"], , drop = FALSE]
dneg <- cd3n[cd3n[[ch["CD19"]]] <= th["CD19"] & cd3n[[ch["CD14"]]] <= th["CD14"], , drop = FALSE]

# A — Mononuclear (FSC-A vs SSC-A); box around the FSC-gated cluster
b <- base_plot(e, "FSC-A", "SSC-A")
bx <- bbox(e, "FSC-A", "SSC-A", e[["FSC-A"]] > FSC_MIN, b$xl, b$yl)
p <- LTL(GATE(b$p, bx), b$xl, b$yl, paste0("Mononuclear cells  ", pc1(nrow(mono), nrow(e))))
sv(p, "A.png")

# B — CD45+ (FSC-A x, CD45 y)
b <- base_plot(mono, "FSC-A", ch["CD45"])
bx <- bbox(mono, "FSC-A", ch["CD45"], mono[[ch["CD45"]]] > th["CD45"], b$xl, b$yl)
p <- LTL(GATE(b$p, bx), b$xl, b$yl, paste0("CD45+  ", pc1(nrow(cd45), nrow(mono))))
sv(p, "B.png")

# C — CD3+ / CD3- (FSC-A x, CD3 y)
b <- base_plot(cd45, "FSC-A", ch["CD3"])
bp <- bbox(cd45, "FSC-A", ch["CD3"], cd45[[ch["CD3"]]] > th["CD3"], b$xl, b$yl)
bn <- bbox(cd45, "FSC-A", ch["CD3"], cd45[[ch["CD3"]]] <= th["CD3"], b$xl, b$yl)
p <- GATE(GATE(b$p, bp), bn)
p <- LTL(p, b$xl, b$yl, paste0("CD3+  ", pc1(nrow(cd3p), nrow(cd45))))
p <- LBR(p, b$xl, b$yl, paste0("CD3-  ", pc1(nrow(cd3n), nrow(cd45))))
sv(p, "C.png")

# D — CD4+ / CD8+ (CD8 x, CD4 y)
b <- base_plot(cd3p, ch["CD8"], ch["CD4"])
m4 <- cd3p[[ch["CD4"]]] > th["CD4"] & cd3p[[ch["CD8"]]] <= th["CD8"]
m8 <- cd3p[[ch["CD8"]]] > th["CD8"] & cd3p[[ch["CD4"]]] <= th["CD4"]
b4 <- bbox(cd3p, ch["CD8"], ch["CD4"], m4, b$xl, b$yl); b8 <- bbox(cd3p, ch["CD8"], ch["CD4"], m8, b$xl, b$yl)
p <- GATE(GATE(b$p, b4), b8)
p <- LTL(p, b$xl, b$yl, paste0("CD4+  ", pc1(sum(m4), nrow(cd3p))))
p <- LBR(p, b$xl, b$yl, paste0("CD8+  ", pc1(sum(m8), nrow(cd3p))))
sv(p, "D.png")

# E — CD14+ / CD19-CD14- (CD19 x, CD14 y)
b <- base_plot(cd3n, ch["CD19"], ch["CD14"])
m14 <- cd3n[[ch["CD14"]]] > th["CD14"]
md  <- cd3n[[ch["CD19"]]] <= th["CD19"] & cd3n[[ch["CD14"]]] <= th["CD14"]
b14 <- bbox(cd3n, ch["CD19"], ch["CD14"], m14, b$xl, b$yl); bd <- bbox(cd3n, ch["CD19"], ch["CD14"], md, b$xl, b$yl)
p <- GATE(GATE(b$p, b14), bd)
p <- LTL(p, b$xl, b$yl, paste0("CD14+  ", pc1(nrow(cd14), nrow(cd3n))))
p <- LBR(p, b$xl, b$yl, paste0("CD19- CD14-  ", pc1(nrow(dneg), nrow(cd3n))))
sv(p, "E.png")

# F — NK CD16+ (CD19 x, CD16 y)
b <- base_plot(dneg, ch["CD19"], ch["CD16"])
m16 <- dneg[[ch["CD16"]]] > th["CD16"] & dneg[[ch["CD19"]]] <= th["CD19"]
b16 <- bbox(dneg, ch["CD19"], ch["CD16"], m16, b$xl, b$yl)
p <- LTL(GATE(b$p, b16), b$xl, b$yl, paste0("NK CD16+  ", pc1(sum(m16), nrow(dneg))))
sv(p, "F.png")

# G — CD64+CD11b+ (CD64 x, CD11b y)
b <- base_plot(cd14, ch["CD64"], ch["CD11b"])
mg <- cd14[[ch["CD64"]]] > th["CD64"] & cd14[[ch["CD11b"]]] > th["CD11b"]
bg <- bbox(cd14, ch["CD64"], ch["CD11b"], mg, b$xl, b$yl)
p <- LTL(GATE(b$p, bg), b$xl, b$yl, paste0("CD64+ CD11b+  ", pc1(sum(mg), nrow(cd14))))
sv(p, "G.png")

message("OK gating_dots A-G (Inflam_Noact_3, boxed)")
