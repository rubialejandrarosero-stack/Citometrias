# R/42_fig_cd45_lines.R — live-cell composition over time: CD45- (spheroid) vs CD45+ (PBMC),
# as % of live cells, one line each in the same panel, faceted by condition. Line = mean of the
# 2 donors, points = donors. Replaces the viability bars with a spheroid-vs-immune comparison.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(scales)})
OUT <- "output/figuras"
POP_COL <- c("CD45- (spheroid)" = "#2c7fb8", "CD45+ (PBMC)" = "#e6550d")

d <- read_csv("output/tidy/cd45_composition.csv", show_col_types = FALSE) %>%
  mutate(cond = factor(paste0(environment, " · ", condition),
                       levels = c("Basal · Resting", "Basal · Activated",
                                  "Inflammatory · Resting", "Inflammatory · Activated")),
         population = factor(population, names(POP_COL)))
s <- d %>% group_by(cond, time, population) %>% summarise(M = mean(value), .groups = "drop")

p <- ggplot(s, aes(time, M, colour = population)) +
  geom_line(linewidth = 1.1) +
  geom_point(data = d, aes(time, value, colour = population), size = 1.6, alpha = 0.7,
             position = position_dodge(width = 2)) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = POP_COL) +
  scale_x_continuous(breaks = c(24, 48, 96)) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ cond, nrow = 1) +
  labs(title = "Live-cell composition: spheroid (CD45-) vs immune (CD45+)",
       x = "Time (h)", y = "% of live cells", colour = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
        strip.background = element_rect(fill = "grey92", colour = "grey70"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(fill = NA, colour = "grey75", linewidth = 0.6),
        panel.spacing = unit(1.1, "lines"), legend.position = "top",
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4))
ggsave(file.path(OUT, "Fig_cd45_lines.png"), p, width = 12, height = 4.6, dpi = 300)
ggsave(file.path(OUT, "Fig_cd45_lines.pdf"), p, width = 12, height = 4.6)
message("OK Fig_cd45_lines.png")
