# R/22_fig_citometria_composicion.R — multi-panel composition of infiltrating immune
# populations (absolute counts) over co-culture days, 4 conditions distinguished by
# colour (basal = black, inflammatory = red) and symbol (resting = open diamond,
# activated = filled square). Style follows the reference panel B.
source("R/fig_theme.R")

grp_lv <- c("Basal, resting", "Basal, activated", "Inflammatory, resting", "Inflammatory, activated")
grp_col <- c("Basal, resting" = "black", "Basal, activated" = "black",
             "Inflammatory, resting" = "#d62728", "Inflammatory, activated" = "#d62728")
grp_shape <- c("Basal, resting" = 5, "Basal, activated" = 15,
               "Inflammatory, resting" = 5, "Inflammatory, activated" = 15)

pops <- c(n_CD45 = "Leukocytes (CD45+)", n_CD3 = "Lymphocytes CD3+", n_CD4 = "Lymphocytes CD4+",
          n_Mono_CD14 = "Monocytes/Macrophages (CD14+)", n_CD8 = "Lymphocytes CD8+", n_NK_CD16 = "NK cells")

d <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(day = recode(Tiempo, `24` = 1, `48` = 2, `96` = 4),
         environment = ifelse(Condicion == "INF", "Inflammatory", "Basal"),
         activation = ifelse(Activacion == "ACT", "activated", "resting"),
         group = factor(paste0(environment, ", ", activation), levels = grp_lv)) %>%
  select(day, group, all_of(names(pops))) %>%
  pivot_longer(all_of(names(pops)), names_to = "pop", values_to = "count") %>%
  mutate(pop = factor(recode(pop, !!!pops), levels = unname(pops)))

s <- d %>% filter(!is.na(count)) %>%
  group_by(pop, group, day) %>%
  summarise(M = mean(count), SE = sd(count) / sqrt(n()), .groups = "drop")

p <- ggplot(s, aes(day, M, colour = group, shape = group, group = group)) +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.12, linewidth = 0.5) +
  geom_point(size = 2.6, stroke = 1, fill = "white") +
  facet_wrap(~pop, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = grp_col) +
  scale_shape_manual(values = grp_shape) +
  scale_x_continuous(breaks = c(1, 2, 4)) +
  labs(title = "Infiltrating immune populations in tumour spheroids",
       subtitle = "Absolute cell counts (mean ± SEM, n = 4 donors)",
       x = "Co-culture time (days)", y = "Cell count (absolute)", colour = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(face = "bold", size = 13),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35"),
        legend.position = "bottom", legend.text = element_text(size = 11))

save_fig(p, "Fig6_composition_counts", w = 9, h = 11)
message("OK Fig6_composition_counts")
