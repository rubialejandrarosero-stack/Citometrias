# R/28_fig_citometria_unificada_portiempo.R — unified dot plot of CD45+ subpopulations
# (CD4, CD8, CD14, NK) faceted by co-culture time (24/48/96 h) in a single image.
# X = 4 conditions (Basal resting/activated | Inflammatory resting/activated), populations
# dodged by colour, individual donors + mean ± SD, on PERCENTAGES (% of CD45+).
# Vertical legend on the right.
source("R/fig_theme.R")
suppressPackageStartupMessages(library(scales))

pop_lv  <- c("CD4+ T cells", "CD8+ T cells", "Macrophages (CD14+)", "NK cells (CD16+)")
pop_pal <- c("CD4+ T cells" = "#2171b5", "CD8+ T cells" = "#31a354",
             "Macrophages (CD14+)" = "#e6550d", "NK cells (CD16+)" = "#756bb1")
cols <- c(pctCD45_CD4 = "CD4+ T cells", pctCD45_CD8 = "CD8+ T cells",
          pctCD45_Mono = "Macrophages (CD14+)", pctCD45_NK = "NK cells (CD16+)")

cx <- c("Basal.Resting" = 1, "Basal.Activated" = 2,
        "Inflammatory.Resting" = 4, "Inflammatory.Activated" = 5)
pop_off <- setNames(c(-0.27, -0.09, 0.09, 0.27), pop_lv)
cond_short <- c("Resting" = "PBMC resting", "Activated" = "PBMC activated")

d <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(environment = ifelse(Condicion == "INF", "Inflammatory", "Basal"),
         cond = ifelse(Activacion == "ACT", "Activated", "Resting"),
         timef = factor(paste0(Tiempo, " h"), levels = c("24 h", "48 h", "96 h"))) %>%
  pivot_longer(all_of(names(cols)), names_to = "pop", values_to = "value") %>%
  mutate(pop = factor(recode(pop, !!!cols), levels = pop_lv),
         base = cx[paste0(environment, ".", cond)], xpos = base + pop_off[as.character(pop)])
s <- d %>% group_by(timef, environment, cond, pop) %>%
  summarise(M = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(base = cx[paste0(environment, ".", cond)], xpos = base + pop_off[as.character(pop)])

# place env labels above the tallest of (points, mean + SD) per timepoint
env_df <- s %>% group_by(timef) %>% summarise(ysd = max(M + SD, na.rm = TRUE), .groups = "drop") %>%
  left_join(d %>% group_by(timef) %>% summarise(yp = max(value, na.rm = TRUE), .groups = "drop"),
            by = "timef") %>%
  mutate(ym = pmax(ysd, yp)) %>% ungroup() %>%
  tidyr::crossing(data.frame(x = c(1.5, 4.5), lab = c("Basal", "Inflammatory"))) %>%
  mutate(y = ym * 1.10)

p <- ggplot() +
  geom_point(data = d, aes(xpos, value, colour = pop), shape = 1, size = 1.8, stroke = 0.9,
             position = position_jitter(width = 0.04, height = 0)) +
  geom_errorbar(data = s, aes(xpos, ymin = M - SD, ymax = M + SD, colour = pop), width = 0.1, linewidth = 0.6) +
  geom_errorbar(data = s, aes(xpos, ymin = M, ymax = M, colour = pop), width = 0.16, linewidth = 0.9) +
  geom_vline(xintercept = 3, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_text(data = env_df, aes(x, y, label = lab), inherit.aes = FALSE, fontface = "bold", size = 5.5) +
  facet_wrap(~timef, nrow = 1, scales = "free") +
  scale_colour_manual(values = pop_pal) +
  scale_x_continuous(breaks = cx, labels = cond_short[sub(".*\\.", "", names(cx))]) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = NULL, y = "% of CD45+ leukocytes", colour = NULL) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme_classic(base_size = 17) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.title.y = element_text(face = "bold", size = 18),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(face = "bold", size = 16),
        legend.position = "right", legend.text = element_text(size = 13),
        panel.grid = element_blank())

save_fig(p, "Fig6_unified_por_tiempo", w = 17, h = 6.5)
message("OK Fig6_unified_por_tiempo")
