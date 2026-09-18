# R/13_fig_citometria.R — Figure 6: infiltrating immune subpopulations (flow cytometry)
source("R/fig_theme.R")

d <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(environment = ifelse(Condicion == "INF", "Inflammatory", "Basal"),
         pbmc = ifelse(Activacion == "ACT", "Activated", "Resting"),
         group = make_group(environment, pbmc),
         time = Tiempo)

plot_pops <- function(cols, labs_map, title, subt, ylab, file, ncol = 3, w = 11, h = 6) {
  long <- d %>% select(group, time, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "pop", values_to = "val") %>%
    mutate(pop = factor(recode(pop, !!!labs_map), levels = unname(labs_map)))
  s <- long %>% filter(!is.na(val)) %>%
    group_by(pop, group, time) %>% summarise(M = mean(val), SE = sem(val), .groups = "drop")
  p <- ggplot(s, aes(time, M, colour = group, group = group)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 3, linewidth = 0.45) +
    facet_wrap(~pop, scales = "free_y", ncol = ncol) +
    scale_colour_manual(values = cond_pal, drop = TRUE) +
    scale_x_continuous(breaks = c(24, 48, 96)) +
    labs(title = title, subtitle = subt, x = "Co-culture time (h)", y = ylab) +
    guides(colour = guide_legend(nrow = 1)) +
    theme_frontiers() +
    theme(legend.text = element_text(margin = margin(r = 14, l = 3)))
  save_fig(p, file, w = w, h = h)
  message("OK ", file)
}

# --- 6A: main immune subpopulations (% of CD45+) ---
plot_pops(
  cols = c("pctCD45_CD3","pctCD45_CD4","pctCD45_CD8","pctCD45_Mono","pctCD45_B","pctCD45_NK"),
  labs_map = c(pctCD45_CD3="CD3+ T cells", pctCD45_CD4="CD4+ T cells", pctCD45_CD8="CD8+ T cells",
               pctCD45_Mono="Monocytes (CD14+)", pctCD45_B="B cells (CD19+)", pctCD45_NK="NK cells (CD16+)"),
  title = "Infiltrating immune subpopulations within tumor spheroids",
  subt = "Flow cytometry, % of CD45+ leukocytes (mean ± SEM, n = 4 donors)",
  ylab = "% of CD45+ leukocytes", file = "Fig6A_subpopulations")

# --- 6B: functional activation (% of parent) ---
plot_pops(
  cols = c("pctPar_CD4HLADR_de_CD4","pctPar_CD8HLADR_de_CD8",
           "pctPar_CD64HLADR_de_Mono","pctPar_CD64CD11b_de_Mono"),
  labs_map = c(pctPar_CD4HLADR_de_CD4="CD4+ HLA-DR+ (% of CD4+)",
               pctPar_CD8HLADR_de_CD8="CD8+ HLA-DR+ (% of CD8+)",
               pctPar_CD64HLADR_de_Mono="CD64+ HLA-DR+ (% of CD14+)",
               pctPar_CD64CD11b_de_Mono="CD64+ CD11b+ (% of CD14+)"),
  title = "Functional activation of infiltrating immune cells",
  subt = "Flow cytometry, activation markers (mean ± SEM, n = 4 donors)",
  ylab = "% of parent population", file = "Fig6B_activation", ncol = 4, w = 12, h = 3.8)
