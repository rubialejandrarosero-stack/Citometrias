# R/34_fig_activacion_dotplot.R — CD64+CD11b+ macrophage activation as quantification
# dot plots (individual donors + mean ± SD + significance), ONE figure per environment
# (Basal / Inflammatory), to complement the representative biaxial figure. New palette
# (teal = resting, purple = activated) to distinguish from the infiltration figures.
source("R/fig_theme.R")
suppressPackageStartupMessages({library(lmerTest); library(emmeans); library(scales)})

POP <- "pctPar_CD64CD11b_de_Mono"
cond_pal <- c("Resting" = "#1B9E77", "Activated" = "#7570B3")   # teal / purple (Dark2)
cond_off <- c("Resting" = -0.19, "Activated" = 0.19)
cond_short <- c("Resting" = "PBMC resting", "Activated" = "PBMC activated")
stars <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", NA))))
tidx <- function(t) match(t, c(24, 48, 96))

theme_prism <- function(base = 18) theme_classic(base_size = base) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = base - 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base - 4),
        axis.title.y = element_text(face = "bold", size = base + 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base + 4),
        legend.position = "none", panel.grid = element_blank())

d0 <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(environment = factor(ifelse(Condicion == "INF", "Inflammatory", "Basal"),
                              levels = c("Basal", "Inflammatory")),
         cond = factor(ifelse(Activacion == "ACT", "Activated", "Resting"), levels = c("Resting", "Activated")),
         timef = factor(Tiempo), donor = factor(Donante), time = Tiempo, value = .data[[POP]]) %>%
  filter(!is.na(value))

m <- lmer(value ~ environment * cond * timef + (1 | donor), data = d0)
ph <- as.data.frame(pairs(emmeans(m, ~ cond | environment * timef), adjust = "tukey")) %>%
  mutate(environment = as.character(environment), time = as.numeric(as.character(timef)), lab = stars(p.value))

fig_env <- function(env) {
  d <- d0 %>% filter(environment == env) %>% mutate(xpos = tidx(time) + cond_off[as.character(cond)])
  s <- d %>% group_by(time, cond) %>% summarise(M = mean(value), SD = sd(value), .groups = "drop") %>%
    mutate(xpos = tidx(time) + cond_off[as.character(cond)])
  gmax <- d %>% group_by(time) %>% summarise(g = max(value), .groups = "drop")
  br <- ph %>% filter(environment == env, !is.na(lab)) %>% left_join(gmax, by = "time") %>%
    mutate(x = tidx(time) + cond_off[["Resting"]], xend = tidx(time) + cond_off[["Activated"]], y = g * 1.05)
  pmax <- max(d$value, na.rm = TRUE); maxc <- if (nrow(br)) max(br$y, pmax) else pmax
  time_y <- maxc * 1.09; ytop <- time_y * 1.06
  axis_df <- s %>% distinct(xpos, cond) %>% mutate(lab = cond_short[as.character(cond)]) %>% arrange(xpos)
  time_df <- distinct(d, time) %>% mutate(x = tidx(time), y = time_y)

  p <- ggplot() +
    geom_point(data = d, aes(xpos, value, colour = cond), shape = 1, size = 2.6, stroke = 1.1,
               position = position_jitter(width = 0.05, height = 0)) +
    geom_errorbar(data = s, aes(xpos, ymin = M - SD, ymax = M + SD, colour = cond), width = 0.13, linewidth = 0.8) +
    geom_errorbar(data = s, aes(xpos, ymin = M, ymax = M, colour = cond), width = 0.26, linewidth = 1.1) +
    geom_text(data = time_df, aes(x, y, label = paste0(time, " h")), inherit.aes = FALSE, fontface = "bold", size = 5.4) +
    scale_colour_manual(values = cond_pal) +
    scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.02, 0.02)), limits = c(0, ytop)) +
    labs(title = env, x = NULL, y = "CD64+CD11b+ (% of CD14+)") +
    theme_prism()
  if (nrow(br)) {
    p <- p + geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.5) +
      geom_text(data = br, aes(x = (x + xend) / 2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.1, size = 6)
  }
  save_fig(p, paste0("Fig6act_CD64CD11b_", tolower(env)), w = 6.5, h = 6)
  message("OK ", env)
}
fig_env("Basal")
fig_env("Inflammatory")
