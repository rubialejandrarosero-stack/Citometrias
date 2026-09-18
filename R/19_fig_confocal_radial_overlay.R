# R/19_fig_confocal_radial_overlay.R — overlaid radial profiles (one curve per
# population) + mean radial position marked; plus a summary table of mean positions.
source("R/fig_theme.R")

pal <- c("PBMC (cell tracker blue)" = "#1f8fff", "A549 (RFP)" = "#ff2b2b",
         "MRC-5 (GFP)" = "#18d63a")

raw <- read_csv("output/tidy/confocal_radial.csv", show_col_types = FALSE)

d <- raw %>%
  mutate(cell = factor(recode(cell_type, MRC5 = "MRC-5 (GFP)", A549 = "A549 (RFP)",
                              PBMC = "PBMC (cell tracker blue)"),
                       levels = c("MRC-5 (GFP)", "A549 (RFP)", "PBMC (cell tracker blue)")),
         condition = factor(recode(condition, Control = "Control", Resting = "Resting PBMC",
                                   Activated = "Activated PBMC"),
                            levels = c("Control", "Resting PBMC", "Activated PBMC")),
         environment = factor(environment, levels = c("Basal", "Inflammatory"),
                              labels = c("Basal", "Inflammatory")))

# mean radial position (weighted by cell count) per population/condition
means <- d %>% group_by(time, environment, condition, cell) %>%
  summarise(mean_pos = sum(dist_pct * count) / sum(count), .groups = "drop")

# ---- summary table (all timepoints) ----
tbl <- means %>%
  mutate(mean_pos = round(mean_pos, 1),
         cell = recode(as.character(cell), "MRC-5 (GFP)" = "MRC5",
                       "A549 (RFP)" = "A549", "PBMC (cell tracker blue)" = "PBMC")) %>%
  tidyr::pivot_wider(names_from = cell, values_from = mean_pos) %>%
  arrange(time, environment, condition)
dir.create("output/estadistica", showWarnings = FALSE, recursive = TRUE)
write_csv(tbl, "output/estadistica/confocal_radial_mean_position.csv")
message("Tabla: output/estadistica/confocal_radial_mean_position.csv")

theme_radial <- theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.7, fill = NA),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 13),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.title = element_text(face = "bold", size = 17),
        axis.title.x = element_text(face = "bold", size = 21),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(colour = "grey35", size = 12),
        legend.position = "bottom", legend.text = element_text(size = 14))

for (tt in c(24, 48, 96)) {
  p <- ggplot(filter(d, time == tt), aes(dist_pct, count, colour = cell, fill = cell)) +
    geom_area(position = "identity", alpha = 0.28, linewidth = 0) +
    geom_line(linewidth = 1) +
    geom_vline(data = filter(means, time == tt),
               aes(xintercept = mean_pos, colour = cell), linetype = "dashed",
               linewidth = 0.7, show.legend = FALSE) +
    geom_vline(xintercept = 100, linetype = "dotted", colour = "grey45", linewidth = 0.5) +
    facet_grid(condition ~ environment) +
    scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = seq(0, 120, 40)) +
    coord_cartesian(xlim = c(0, 120)) +
    labs(title = paste0("Radial profiles of cell populations (", tt, " h)"),
         subtitle = "Coloured dashed line = mean radial position; grey dotted = spheroid edge. n = 1 (descriptive)",
         x = "Distance from spheroid center (% of radius)", y = "Number of cells",
         colour = NULL, fill = NULL) +
    theme_radial
  # bare version for the vertical montage: title = just the time; X-axis title only
  # on the 96 h block (shown once at the very bottom).
  xt <- if (tt == 96) "Distance from spheroid center (% of radius)" else NULL
  pbare <- p + labs(title = paste0(tt, " h"), subtitle = NULL, x = xt) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0, size = 24))
  save_fig(pbare, paste0("Fig4_radial_bare_", tt, "h"), w = 9.5, h = if (tt == 96) 7.7 else 7.2)
  # clean body version (single timepoint): title = just the time, no subtitle, legend kept
  save_fig(p + labs(title = paste0(tt, " h"), subtitle = NULL),
           paste0("Fig4_radial_body_", tt, "h"), w = 9.5, h = 8.2)
  message("OK Fig4_radial_overlay_", tt, "h")
}

# --- standalone legend (vertical) for the montage ---
suppressPackageStartupMessages({library(grid); library(gridExtra)})
pl <- ggplot(filter(d, time == 24), aes(dist_pct, count, colour = cell, fill = cell)) +
  geom_area(position = "identity", alpha = 0.28) + geom_line(linewidth = 1) +
  scale_fill_manual(values = pal) + scale_colour_manual(values = pal) +
  guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1)) +
  labs(fill = NULL, colour = NULL) +
  theme(legend.text = element_text(size = 17), legend.title = element_blank(),
        legend.key.size = unit(1.0, "cm"),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        legend.margin = margin(8, 12, 8, 8))
gt <- ggplotGrob(pl)
lg <- gt$grobs[[which(vapply(gt$grobs, function(g) g$name, "") == "guide-box")]]
png("output/figuras/Fig4_radial_legend.png", width = 720, height = 420, res = 150, bg = "white")
grid.newpage(); grid.draw(lg); dev.off()
message("OK Fig4_radial_legend")
