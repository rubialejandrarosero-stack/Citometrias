# R/30_dotplots_representativos.R — representative biaxial density dot plots (pseudocolor)
# with gate + % positive, for the manuscript's activation/gating panel.
# Grid: 4 conditions (environment x immune response) x 4 activation biaxials, at 96 h.
# Representative donor per condition = closest to the condition mean of the 4 activation %.
suppressPackageStartupMessages({library(flowCore); library(ggplot2); library(dplyr); library(readr)})

DIR <- "data/raw/Infiltracion"; REF <- file.path(DIR, "Reference Group")
OUT <- "output/figuras/dotplots"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ch <- c(CD45="APC-Cy7-A", CD3="PE-Cy7-A", CD4="BV421-A", CD8="BV510-A", CD14="PE-A",
        CD16="FITC-A", CD19="StarBrightViolet 475-A", CD64="PerCP-Cy5.5-A",
        CD11b="APC-A", HLADR="Pacific Blue-A")
FSC_MIN <- 5e5
rd <- function(f) read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
lg <- estimateLogicle(rd(file.path(DIR, "48h/INF_ACT/donante_3.fcs")), channels = unname(ch))
tdf <- function(f) { e <- as.data.frame(exprs(transform(rd(f), lg))); e[e[["FSC-A"]] > FSC_MIN, , drop = FALSE] }
un <- tdf(file.path(REF, "Unstained (Cells).fcs"))
th <- setNames(sapply(ch, function(c) as.numeric(quantile(un[[c]], 0.999, na.rm = TRUE))), names(ch))
g <- function(d, m, op = `>`) d[op(d[[ch[[m]]]], th[m]), , drop = FALSE]
flow_ramp <- colorRampPalette(c("#000091", "#0000FF", "#00C0FF", "#00FF60",
                                "#B0FF00", "#FFD000", "#FF3000", "#B00000"))

# biaxial density plot; gate = upper-right (mx+ my+); pct = gate / denom
biax <- function(e, mx, my, denom, file) {
  d <- data.frame(x = e[[ch[[mx]]]], y = e[[ch[[my]]]])
  d$col <- densCols(d$x, d$y, colramp = flow_ramp, nbin = 256)
  ingate <- sum(d$x > th[mx] & d$y > th[my]); pct <- round(ingate / denom * 100, 1)
  qx <- quantile(d$x, c(0.002, 0.9997), na.rm = TRUE); qy <- quantile(d$y, c(0.002, 0.9997), na.rm = TRUE)
  xlim <- c(min(qx[1], th[mx] - 0.4), qx[2]); ylim <- c(min(qy[1], th[my] - 0.4), qy[2])
  p <- ggplot(d, aes(x, y)) +
    geom_point(colour = d$col, size = 0.2, shape = 16) +
    geom_vline(xintercept = th[mx], colour = "grey35", linewidth = 0.3) +
    geom_hline(yintercept = th[my], colour = "grey35", linewidth = 0.3) +
    annotate("rect", xmin = th[mx], xmax = xlim[2], ymin = th[my], ymax = ylim[2],
             fill = NA, colour = "black", linewidth = 0.5) +
    annotate("text", x = xlim[2], y = ylim[2], label = paste0(pct, "%"),
             hjust = 1.1, vjust = 1.4, fontface = "bold", size = 5) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = mx, y = gsub("HLADR", "HLA-DR", my)) +
    theme_classic(base_size = 15) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_text(face = "bold", size = 16),
          plot.margin = margin(6, 6, 6, 6),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6))
  ggsave(file.path(OUT, file), p, width = 3.4, height = 3.4, dpi = 200)
}

# --- representative donor per condition (96 h) ---
csv <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  filter(Tiempo == 96)
actcols <- c("pctPar_CD4HLADR_de_CD4","pctPar_CD8HLADR_de_CD8",
             "pctPar_CD64HLADR_de_Mono","pctPar_CD64CD11b_de_Mono")
pick_rep <- function(cond, act) {
  sub <- csv %>% filter(Condicion == cond, Activacion == act)
  M <- as.matrix(sub[, actcols]); z <- scale(M)
  z[is.na(z)] <- 0; ctr <- colMeans(z)
  sub$Donante[which.min(rowSums((z - matrix(ctr, nrow(z), ncol(z), byrow = TRUE))^2))]
}

conds <- list(
  list(key="Bas_rest", folder="NO_INF_NO_ACT", cond="NO_INF", act="NO_ACT"),
  list(key="Bas_act",  folder="NO_INF_ACT",    cond="NO_INF", act="ACT"),
  list(key="Inf_rest", folder="INF_NO_ACT",    cond="INF",    act="NO_ACT"),
  list(key="Inf_act",  folder="INF_ACT",       cond="INF",    act="ACT"))

# Only CD64/CD11b (macrophage activation) is shown: HLA-DR-based gates were removed
# because HLA-DR does not resolve in this violet-heavy spectral panel (see R/33).
for (cc in conds) {
  don <- pick_rep(cc$cond, cc$act)
  e <- tdf(file.path(DIR, "96h", cc$folder, paste0("donante_", don, ".fcs")))
  cd45 <- g(e, "CD45")
  cd3n <- cd45[cd45[[ch[["CD3"]]]] <= th["CD3"], , drop = FALSE]; cd14 <- g(cd3n, "CD14")
  biax(cd14, "CD64", "CD11b", nrow(cd14), paste0(cc$key, "_CD64CD11b.png"))
  message("OK ", cc$key, " (donante ", don, ")")
}
