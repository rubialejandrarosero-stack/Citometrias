# R/fig_theme.R — unified figure style for the manuscript (Frontiers)
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(readr)})

# Consistent condition scheme across ALL figures:
#   environment: colour family (Basal = blue, Inflammatory = red), controls = grey
#   PBMC state:  shade (resting = light, activated = dark)
cond_levels <- c("Control (basal)","Basal + resting PBMC","Basal + activated PBMC",
                 "Control (inflammatory)","Inflammatory + resting PBMC","Inflammatory + activated PBMC")
cond_pal <- c(
  "Control (basal)"                = "#bdbdbd",
  "Basal + resting PBMC"           = "#9ecae1",
  "Basal + activated PBMC"         = "#2171b5",
  "Control (inflammatory)"         = "#737373",
  "Inflammatory + resting PBMC"    = "#fdae6b",
  "Inflammatory + activated PBMC"  = "#cb181d")

make_group <- function(environment, pbmc) {
  env <- ifelse(environment == "Inflammatory", "Inflammatory", "Basal")
  lab <- dplyr::case_when(
    pbmc == "None"      ~ paste0("Control (", tolower(env), ")"),
    pbmc == "Resting"   ~ paste0(env, " + resting PBMC"),
    pbmc == "Activated" ~ paste0(env, " + activated PBMC"))
  factor(lab, levels = cond_levels)
}

theme_frontiers <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
          plot.title = element_text(face = "bold", size = base + 2, hjust = 0),
          plot.subtitle = element_text(colour = "grey35", size = base - 1),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom", legend.title = element_blank(),
          legend.key.size = unit(0.9, "lines"), legend.text = element_text(size = base - 1))
}

sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# save both PNG (300 dpi) and vector PDF
save_fig <- function(p, name, w = 8, h = 5) {
  dir.create("output/figuras", recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path("output/figuras", paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path("output/figuras", paste0(name, ".pdf")), p, width = w, height = h)
}
