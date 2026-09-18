# R/33_hladr_resolution.R — does HLA-DR resolve in this panel? In resting samples,
# compare HLA-DR on CD4+ T cells (expected HLA-DR-negative) vs CD14+ monocytes
# (bona-fide HLA-DR-high APCs). If monocytes >> T cells -> marker resolves (any high
# T-cell signal is a threshold problem). If they overlap high -> spreading/technical.
suppressPackageStartupMessages({library(flowCore); library(ggplot2); library(dplyr)})

DIR <- "data/raw/Infiltracion"; REF <- file.path(DIR, "Reference Group")
OUT <- "output/figuras/diagnostics"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ch <- c(CD45="APC-Cy7-A", CD3="PE-Cy7-A", CD4="BV421-A", CD14="PE-A", HLADR="Pacific Blue-A")
FSC_MIN <- 5e5
rd <- function(f) read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
lg <- estimateLogicle(rd(file.path(DIR, "48h/INF_ACT/donante_3.fcs")), channels = unname(ch))
tdf <- function(f) { e <- as.data.frame(exprs(transform(rd(f), lg))); e[e[["FSC-A"]] > FSC_MIN, , drop = FALSE] }
un <- tdf(file.path(REF, "Unstained (Cells).fcs"))
th <- setNames(sapply(ch, function(c) as.numeric(quantile(un[[c]], 0.999, na.rm = TRUE))), names(ch))
g <- function(d, m, op = `>`) d[op(d[[ch[[m]]]], th[m]), , drop = FALSE]

samples <- list(
  `Basal resting (96 h)`        = "96h/NO_INF_NO_ACT/donante_3.fcs",
  `Inflammatory resting (96 h)` = "96h/INF_NO_ACT/donante_1.fcs",
  `Basal resting (48 h)`        = "48h/NO_INF_NO_ACT/donante_1.fcs")

cat(sprintf("HLA-DR threshold (logicle) = %.2f ; unstained median = %.2f\n\n",
            th["HLADR"], median(un[[ch[["HLADR"]]]])))
cat(sprintf("%-28s %-14s %6s %8s %10s\n", "sample", "population", "n", "%HLADR+", "median"))
cat(strrep("-", 72), "\n")
dens <- list()
for (nm in names(samples)) {
  e <- tdf(file.path(DIR, samples[[nm]]))
  cd45 <- g(e, "CD45"); cd3 <- g(cd45, "CD3")
  cd4  <- g(cd3, "CD4")
  cd3n <- cd45[cd45[[ch[["CD3"]]]] <= th["CD3"], , drop = FALSE]; mono <- g(cd3n, "CD14")
  for (p in list(list("CD4+ T cells", cd4), list("CD14+ monocytes", mono))) {
    v <- p[[2]][[ch[["HLADR"]]]]
    cat(sprintf("%-28s %-14s %6d %7.1f%% %10.2f\n", nm, p[[1]], length(v),
                mean(v > th["HLADR"]) * 100, median(v)))
    dens[[length(dens) + 1]] <- data.frame(sample = nm, population = p[[1]], hladr = v)
  }
}
dd <- bind_rows(dens)
dd$population <- factor(dd$population, levels = c("CD4+ T cells", "CD14+ monocytes"))

p <- ggplot(dd, aes(hladr, colour = population, fill = population)) +
  geom_density(alpha = 0.25, linewidth = 0.9) +
  geom_vline(xintercept = th["HLADR"], linetype = "dashed", colour = "grey30") +
  geom_vline(xintercept = median(un[[ch[["HLADR"]]]]), linetype = "dotted", colour = "black") +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("CD4+ T cells" = "#2171b5", "CD14+ monocytes" = "#e6550d")) +
  scale_fill_manual(values = c("CD4+ T cells" = "#2171b5", "CD14+ monocytes" = "#e6550d")) +
  labs(x = "HLA-DR (Pacific Blue, logicle)", y = "Density", colour = NULL, fill = NULL,
       title = "HLA-DR resolution: T cells vs monocytes",
       subtitle = "dashed = p99.9 unstained threshold ; dotted = unstained median") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUT, "hladr_resolution.png"), p, width = 7, height = 7.5, dpi = 150)
message("Guardado: ", file.path(OUT, "hladr_resolution.png"))
