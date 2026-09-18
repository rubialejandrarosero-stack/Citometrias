# R/12_fig_confocal.R — Figure 4: confocal quantification (PBMC infiltration, A549, MRC-5)
# Sample 1-6 per timepoint = the 6 conditions in order (confirmed by author).
source("R/fig_theme.R")

map <- tibble::tibble(
  sample = 1:6,
  environment = c("Basal","Basal","Basal","Inflammatory","Inflammatory","Inflammatory"),
  pbmc = c("None","Resting","Activated","None","Resting","Activated"))

d <- read_csv("output/tidy/confocal_tidy.csv", show_col_types = FALSE) %>%
  left_join(map, by = "sample") %>%
  mutate(group = make_group(environment, pbmc),
         celltype = recode(celltype,
                           PBMC = "PBMC infiltration",
                           A549 = "A549 tumor cells",
                           MRC5 = "MRC-5 fibroblasts"),
         celltype = factor(celltype, levels = c("PBMC infiltration","A549 tumor cells","MRC-5 fibroblasts")))

p <- ggplot(d, aes(time, norm_total, colour = group, group = group)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  facet_wrap(~celltype, nrow = 1) +
  scale_colour_manual(values = cond_pal, drop = FALSE) +
  scale_x_continuous(breaks = c(24, 48, 96)) +
  labs(title = "Immune infiltration and tumor/stromal dynamics (confocal microscopy)",
       subtitle = "Normalized fluorescence intensity per condition (single spheroid per condition)",
       x = "Co-culture time (h)", y = "Normalized intensity (a.u.)") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_frontiers() +
  theme(legend.text = element_text(margin = margin(r = 16, l = 3)),
        legend.key.width = unit(1.5, "lines"))

save_fig(p, "Fig4_confocal", w = 12, h = 4.6)
message("OK Fig4_confocal")
