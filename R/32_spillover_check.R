# R/32_spillover_check.R — diagnose spectral spillover/spreading into the HLA-DR
# (Pacific Blue) channel using the single-stain controls. Violet dyes BV421 (CD4),
# BV510 (CD8) and Pacific Blue (HLA-DR) share the 405 nm laser -> candidate spillover.
# For each single stain: among cells positive in their OWN marker, what % cross the
# HLA-DR threshold (i.e. appear falsely HLA-DR+). Thresholds = p99.9 of Unstained.
suppressPackageStartupMessages({library(flowCore); library(ggplot2); library(dplyr)})

DIR <- "data/raw/Infiltracion"; REF <- file.path(DIR, "Reference Group")
OUT <- "output/figuras/diagnostics"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ch <- c(CD45="APC-Cy7-A", CD4="BV421-A", CD8="BV510-A", CD14="PE-A", CD16="FITC-A",
        HLADR="Pacific Blue-A")
FSC_MIN <- 5e5
rd <- function(f) read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
lg <- estimateLogicle(rd(file.path(DIR, "48h/INF_ACT/donante_3.fcs")), channels = unname(ch))
tdf <- function(f) { e <- as.data.frame(exprs(transform(rd(f), lg))); e[e[["FSC-A"]] > FSC_MIN, , drop = FALSE] }
un <- tdf(file.path(REF, "Unstained (Cells).fcs"))
th <- setNames(sapply(ch, function(c) as.numeric(quantile(un[[c]], 0.999, na.rm = TRUE))), names(ch))

# single stains to test as spillover SOURCES into HLA-DR
stains <- list(
  CD4  = list(file = "CD4 BV421 (Cells).fcs",        own = "CD4",   laser = "violet (BV421)"),
  CD8  = list(file = "CD8 BV510 (Cells).fcs",        own = "CD8",   laser = "violet (BV510)"),
  HLADR= list(file = "HLA-DR Pacific Blue (Cells).fcs", own = "HLADR", laser = "violet (Pac Blue)"),
  CD14 = list(file = "CD14 PE (Cells).fcs",          own = "CD14",  laser = "blue-green (PE)"),
  CD16 = list(file = "CD16 FITC (Cells).fcs",        own = "CD16",  laser = "blue (FITC)"))

cat(sprintf("Unstained HLA-DR threshold (logicle) = %.2f\n\n", th["HLADR"]))
cat(sprintf("%-8s %-18s %8s %10s %12s %12s\n",
            "stain", "laser", "n_own+", "%HLADR+", "med.HLADR+", "med.HLADR-"))
cat(strrep("-", 74), "\n")
rows <- list()
for (nm in names(stains)) {
  s <- stains[[nm]]; e <- tdf(file.path(REF, s$file))
  ownpos <- e[e[[ch[[s$own]]]] > th[s$own], , drop = FALSE]
  if (!nrow(ownpos)) next
  pct_hladr <- mean(ownpos[[ch[["HLADR"]]]] > th["HLADR"]) * 100
  med_pos <- median(ownpos[[ch[["HLADR"]]]])
  med_neg <- median(e[e[[ch[[s$own]]]] <= th[s$own], ch[["HLADR"]]], na.rm = TRUE)
  cat(sprintf("%-8s %-18s %8d %9.1f%% %12.2f %12.2f\n",
              nm, s$laser, nrow(ownpos), pct_hladr, med_pos, med_neg))
  rows[[nm]] <- data.frame(stain = nm, laser = s$laser, n = nrow(ownpos),
                           pct_hladr_pos = round(pct_hladr, 1))
}
cat("\nInterpretation: a violet single stain (CD4/CD8) with a HIGH %HLADR+ and med.HLADR+\n",
    ">> med.HLADR- indicates its dye spills/spreads into the HLA-DR channel (false HLA-DR+).\n",
    "A non-violet stain (PE/FITC) should stay near the unstained level (~0.1%).\n", sep = "")

# --- diagnostic biaxials: own marker (x) vs HLA-DR (y) for the violet stains ---
flow_ramp <- colorRampPalette(c("#000091","#0000FF","#00C0FF","#00FF60","#B0FF00","#FFD000","#FF3000","#B00000"))
biax <- function(nm) {
  s <- stains[[nm]]; e <- tdf(file.path(REF, s$file))
  d <- data.frame(x = e[[ch[[s$own]]]], y = e[[ch[["HLADR"]]]])
  d$col <- densCols(d$x, d$y, colramp = flow_ramp, nbin = 200)
  ggplot(d, aes(x, y)) + geom_point(colour = d$col, size = 0.25) +
    geom_vline(xintercept = th[s$own], colour = "grey30") +
    geom_hline(yintercept = th["HLADR"], colour = "grey30") +
    labs(title = paste0(nm, " single stain"), x = s$own, y = "HLA-DR (Pac Blue)") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black"))
}
for (nm in c("CD4","CD8","HLADR"))
  ggsave(file.path(OUT, paste0("spillover_", nm, ".png")), biax(nm), width = 3.6, height = 3.6, dpi = 150)
message("Biaxiales guardados en ", OUT)
