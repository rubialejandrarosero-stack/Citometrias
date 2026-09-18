# R/17_fig_confocal_composicion.R — Figure 4 (confocal): stacked composition.
# % of total fluorescence intensity per cell type (PBMC / A549 / MRC-5), split into two panels
# (Basal | cytokine cocktail). n = 1 spheroid per condition/time -> descriptive (no error bars).
source("R/fig_theme.R")

cell_pal <- c("MRC-5 (GFP)" = "#18d63a", "A549 (KAT)" = "#ff2b2b", "PBMC (blue)" = "#1f8fff")
day_lab  <- c("24" = "1 day", "48" = "2 days", "96" = "4 days")
cond_off <- c("Control" = -0.25, "Resting" = 0, "Activated" = 0.25)
tidx <- function(t) match(t, c(24, 48, 96))

map <- tibble::tibble(sample = 1:6,
  env  = c("Basal", "Basal", "Basal", "Inflammatory", "Inflammatory", "Inflammatory"),
  cond = c("Control", "Resting", "Activated", "Control", "Resting", "Activated"))

d <- read_csv("output/tidy/confocal_tidy.csv", show_col_types = FALSE) %>%
  select(celltype, time, sample, total) %>%
  group_by(celltype) %>% mutate(bn = total / max(total)) %>% ungroup() %>%
  left_join(map, by = "sample") %>%
  group_by(sample, time) %>% mutate(pct = 100 * bn / sum(bn)) %>% ungroup() %>%
  mutate(celltype = recode(celltype, PBMC = "PBMC (blue)", A549 = "A549 (KAT)",
                           MRC5 = "MRC-5 (GFP)"),
         celltype = factor(celltype, levels = names(cell_pal)),
         cond = factor(cond, levels = c("Control", "Resting", "Activated")),
         env2 = factor(env, levels = c("Basal", "Inflammatory"),
                       labels = c("Basal", "TNF α, IL-α, IL-1 β")),
         xpos = tidx(time) + cond_off[as.character(cond)])

time_y  <- 108
time_df <- d %>% distinct(env2, time) %>%
  mutate(x = tidx(time), lab = day_lab[as.character(time)])
axis_df <- d %>% distinct(xpos, cond) %>%
  mutate(lab = recode(as.character(cond), Control = "Control",
                      Resting = "Resting PBMC", Activated = "Activated PBMC")) %>% arrange(xpos)

p <- ggplot(d, aes(xpos, pct, fill = celltype)) +
  geom_col(width = 0.24, colour = "black", linewidth = 0.25) +
  geom_text(data = time_df, aes(x, time_y, label = lab), inherit.aes = FALSE,
            fontface = "bold", size = 6) +
  facet_wrap(~ env2, nrow = 1) +
  scale_fill_manual(values = cell_pal) +
  scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 116), clip = "off") +
  labs(x = NULL, y = "% of total cells", fill = NULL) +
  theme_classic(base_size = 19) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.text = element_text(colour = "black", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.title = element_text(face = "bold", size = 19),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 19),
        panel.spacing = unit(1.3, "lines"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = 16),
        panel.grid = element_blank(), plot.margin = margin(6, 12, 6, 6))

save_fig(p, "Fig4_confocal_composition", w = 13, h = 6.8)
message("OK Fig4_confocal_composition")
