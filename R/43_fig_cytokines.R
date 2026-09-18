# R/43_fig_cytokines.R — key cytokines from the Luminex secretome (pg/mL): CCL2, CCL3, CCL4,
# CXCL9, IL-8, IFN-g, IL-12/IL-23p40, IL-10, IL-1RA, sCD25/sIL-2Ra, VEGF, HGF, FGF-2, G-CSF,
# GM-CSF. Each cytokine is an INDEPENDENT panel with its own plain title, broken down by its 3
# timepoints (24/48/96 h), Basal | Inflammatory groups (thick divider + in-panel label),
# Control/Resting/Activated as dodged HOLLOW bars (colour = condition, black/blue/red) + donor
# points + SD. Arranged 3 x 5, single shared legend at the bottom.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales); library(grid); library(tidyr)
})
OUT <- "output/figuras"
cond_lv  <- c("Control", "Resting", "Activated")
cond_pal <- c(Control = "black", Resting = "#4292c6", Activated = "#e31a1c")
cond_lab <- c("Control spheroid", "Resting PBMCs", "Activated PBMCs")
cond_off <- c(Control = -0.70, Resting = 0, Activated = 0.70)
env_lab  <- c(Basal = "Basal", Inflammatory = "TNF α, IL-α, IL-1 β")
tidx <- function(t) match(t, c(24, 48, 96)) * 2
xbase <- function(env, t) ifelse(env == "Basal", tidx(t), tidx(t) + 8)
KEY <- c("CCL2" = "CCL2", "CCL3" = "CCL3", "CCL4" = "CCL4", "CXCL-9" = "CXCL9",
         "IL-8" = "IL-8", "IFN-gamma" = "IFN-γ", "IL-12/IL-23p40" = "IL-12/IL-23p40",
         "IL-10" = "IL-10", "IL-1RA" = "IL-1RA", "CD25" = "sCD25 / sIL-2Rα",
         "VEGF" = "VEGF", "HGF" = "HGF", "FGF basic/ FGF2/bFGF" = "FGF-2",
         "G-CSF" = "G-CSF", "GM-CSF" = "GM-CSF")
# title colour by biological family, matching the Luminex heatmap taxonomy
FAM_COL <- c(CCL2 = "#1f6fb2", CCL3 = "#1f6fb2", CCL4 = "#1f6fb2", CXCL9 = "#1f6fb2",
             "IL-8" = "#1f6fb2",                                   # Chemokines - blue
             "IFN-γ" = "#111111", "IL-12/IL-23p40" = "#111111",    # Pro-inflammatory - black
             "IL-10" = "#2e8b3d", "IL-1RA" = "#2e8b3d",            # Th2 / regulatory - green
             "sCD25 / sIL-2Rα" = "#111111",                        # T-cell - black
             VEGF = "#e08214", HGF = "#e08214", "FGF-2" = "#e08214",
             "G-CSF" = "#e08214", "GM-CSF" = "#e08214")            # Growth factors - orange

d0 <- read_csv("output/tidy/luminex_tidy.csv", show_col_types = FALSE) %>%
  filter(cytokine %in% names(KEY)) %>%
  mutate(condition = factor(ifelse(pbmc %in% c("None", NA), "Control", pbmc), cond_lv),
         cyt = factor(KEY[cytokine], levels = unname(KEY)),
         environment = factor(environment, c("Basal", "Inflammatory")), timef = factor(time)) %>%
  group_by(cytokine, environment, condition, time) %>% mutate(rep_id = row_number()) %>% ungroup() %>%
  mutate(donor2 = factor(ifelse(condition == "Control", paste0("C", rep_id), as.character(donor))),
         x = xbase(environment, time) + cond_off[as.character(condition)])
s0 <- d0 %>% group_by(cyt, environment, condition, time, x) %>%
  summarise(M = mean(value), SD = sd(value), .groups = "drop")

panel <- function(cy) {
  dd <- filter(d0, cyt == cy); ss <- filter(s0, cyt == cy)
  bar_top <- max(ss$M + replace_na(ss$SD, 0), dd$value, na.rm = TRUE)
  ytop <- bar_top * 1.18

  # in-panel environment labels ("Basal" / "Inflammatory"), centred over each half, near the top
  env_labs <- tibble(x = c(4, 12), lab = c("Basal", "Inflammatory"), y = ytop * 0.98)

  p <- ggplot() +
    geom_col(data = ss, aes(x, M, colour = condition), fill = NA, width = 0.35, linewidth = 2.6) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + SD, colour = condition),
                  width = 0.19, linewidth = 1.6) +
    geom_point(data = dd, aes(x, value, colour = condition), size = 1.6, alpha = 0.8,
               position = position_jitter(width = 0.05, height = 0)) +
    geom_vline(xintercept = 8, colour = "black", linewidth = 1.3) +      # thick divider Basal | cocktail
    geom_text(data = env_labs, aes(x, y, label = lab), inherit.aes = FALSE,
              size = 14, fontface = "italic", colour = "grey30", vjust = 1) +
    scale_colour_manual(values = cond_pal, labels = cond_lab, name = NULL) +
    scale_x_continuous(breaks = c(2, 4, 6, 10, 12, 14), labels = rep(c("24", "48", "96"), 2),
                       limits = c(1.0, 15.0)) +
    scale_y_continuous(labels = label_comma(), limits = c(0, ytop),
                       expand = expansion(mult = c(0, 0.02))) +      # hard floor at 0 -> bars sit on the axis
    labs(x = "Time (h)", y = "Concentration (pg/mL)", title = as.character(cy)) +
    theme_classic(base_size = 46) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.7),
          axis.ticks.x = element_blank(),
          axis.text = element_text(colour = "black", size = 46),
          axis.text.x = element_text(size = 40, angle = 0, hjust = 0.5),
          axis.title.y = element_text(face = "bold", size = 46),
          axis.title.x = element_text(face = "bold", size = 46),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 48,
                                     colour = FAM_COL[as.character(cy)]),
          legend.position = "none", panel.grid = element_blank())

  p
}

cys <- unname(KEY)
NCOL <- 3
nrows_panels <- ceiling(length(cys) / NCOL)

# single shared legend, extracted from an auxiliary plot that does carry one
leg_plot <- panel(cys[1]) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 40), legend.key.size = unit(2.0, "lines"),
        legend.margin = margin(10, 10, 10, 10))
leg_grob <- ggplotGrob(leg_plot)
legend <- leg_grob$grobs[[which(vapply(leg_grob$grobs, function(g) g$name, "") == "guide-box")]]
leg_h <- unit(4.2, "lines")

png(file.path(OUT, "Fig_cytokines_key.png"), width = 12.5 * NCOL, height = 9.0 * nrows_panels + 0.8,
    units = "in", res = 220, bg = "white")
grid.newpage()
heights_u <- rep(unit.c(unit(1, "null"), unit(2.4, "lines")), nrows_panels)
heights_u <- heights_u[seq_len(length(heights_u) - 1)]
heights_u <- unit.c(leg_h, heights_u)
pushViewport(viewport(layout = grid.layout(length(heights_u), NCOL, heights = heights_u)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:NCOL))
grid.draw(legend)
popViewport()
prt <- function(p, r, c) {
  pushViewport(viewport(layout.pos.row = r, layout.pos.col = c, clip = "off"))  # let below-axis daggers bleed out
  print(p, newpage = FALSE); popViewport()
}
panel_rows <- seq(2, by = 2, length.out = nrows_panels)
for (i in seq_along(cys)) {
  r <- panel_rows[(i - 1) %/% NCOL + 1]
  prt(panel(cys[i]), r, (i - 1) %% NCOL + 1)
}
popViewport(); dev.off()
message("OK Fig_cytokines_key.png")

# Individual per-cytokine files (same panel design/size as one cell of the combined grid)
IND_DIR <- file.path(OUT, "cytokines_individual")
dir.create(IND_DIR, showWarnings = FALSE, recursive = TRUE)
slug <- function(s) {
  s <- gsub("α", "alpha", s); s <- gsub("β", "beta", s); s <- gsub("γ", "gamma", s)
  gsub("[^A-Za-z0-9]+", "", s)
}
for (i in seq_along(cys)) {
  fn <- file.path(IND_DIR, paste0("Fig_cytokine_", slug(cys[i]), ".png"))
  p_ind <- panel(cys[i]) +
    theme(legend.position = "top", legend.direction = "horizontal",
          legend.text = element_text(size = 40), legend.key.size = unit(1.8, "lines"),
          legend.margin = margin(1, 1, 1, 1))
  ggsave(fn, p_ind, width = 17.0, height = 12.0, units = "in", dpi = 300, bg = "white")
  message("OK ", fn)
}
