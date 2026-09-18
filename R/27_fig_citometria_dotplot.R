# R/27_fig_citometria_dotplot.R — scatter dot plot (individual donors) + mean ± SD,
# no bars, significance from the linear mixed model (env x activation x time + (1|donor)).
# Population name goes in the X-axis title (no plot title). Same X layout as the rest.
source("R/fig_theme.R")
suppressPackageStartupMessages({library(lmerTest); library(emmeans); library(scales)})

cond2_pal <- c("Resting" = "#4292c6", "Activated" = "#e6550d")
cond_off  <- c("Resting" = -0.19, "Activated" = 0.19)
cond_short <- c("Resting" = "PBMC resting", "Activated" = "PBMC activated")
stars <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", NA))))
tidx  <- function(t) match(t, c(24, 48, 96))
xbase <- function(env, t) ifelse(env == "Basal", tidx(t), tidx(t) + 4)

theme_prism <- function(base = 18) theme_classic(base_size = base) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = base - 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base - 5),
        axis.title.x = element_text(face = "bold", size = base + 4, margin = margin(t = 8)),
        axis.title.y = element_text(face = "bold", size = base + 1),
        legend.position = "none", panel.grid = element_blank())

d0 <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(environment = factor(ifelse(Condicion == "INF", "Inflammatory", "Basal"),
                              levels = c("Basal", "Inflammatory")),
         cond = factor(ifelse(Activacion == "ACT", "Activated", "Resting"),
                       levels = c("Resting", "Activated")),
         timef = factor(Tiempo), donor = factor(Donante), time = Tiempo)

fig_dot <- function(pop, ylab, file) {
  d <- d0 %>% mutate(value = .data[[pop]]) %>% filter(!is.na(value)) %>%
    mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])
  s <- d %>% group_by(environment, time, cond) %>%
    summarise(M = mean(value), SD = sd(value), .groups = "drop") %>%
    mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])

  ph <- tryCatch({
    m <- lmer(value ~ environment * cond * timef + (1 | donor), data = d)
    as.data.frame(pairs(emmeans(m, ~ cond | environment * timef), adjust = "tukey"))
  }, error = function(e) NULL)

  gmax <- d %>% group_by(environment, time) %>%
    summarise(g = max(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(environment = as.character(environment))
  br <- NULL
  if (!is.null(ph)) {
    br <- ph %>% mutate(environment = as.character(environment),
                        time = as.numeric(as.character(timef)), lab = stars(p.value)) %>%
      filter(!is.na(lab)) %>% left_join(gmax, by = c("environment", "time")) %>%
      mutate(x = xbase(environment, time) + cond_off[["Resting"]],
             xend = xbase(environment, time) + cond_off[["Activated"]], y = g * 1.07)
  }
  pmax <- max(d$value, na.rm = TRUE)
  maxc <- if (!is.null(br) && nrow(br)) max(br$y, pmax) else pmax
  time_y <- maxc * 1.10; env_y <- time_y * 1.07; ymax <- env_y * 1.05

  axis_df <- s %>% distinct(xpos, cond) %>% mutate(lab = cond_short[as.character(cond)]) %>% arrange(xpos)
  time_df <- s %>% distinct(environment, time) %>% mutate(x = xbase(environment, time))

  p <- ggplot() +
    geom_point(data = d, aes(xpos, value, colour = cond), shape = 1, size = 2.4, stroke = 1,
               position = position_jitter(width = 0.05, height = 0)) +
    geom_errorbar(data = s, aes(xpos, ymin = M - SD, ymax = M + SD, colour = cond),
                  width = 0.14, linewidth = 0.7) +
    geom_errorbar(data = s, aes(xpos, ymin = M, ymax = M, colour = cond),
                  width = 0.28, linewidth = 1) +
    geom_vline(xintercept = 4, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text(data = time_df, aes(x, time_y, label = paste0(time, " h")),
              inherit.aes = FALSE, fontface = "bold", size = 5) +
    annotate("text", x = 2, y = env_y, label = "Basal", fontface = "bold", size = 7) +
    annotate("text", x = 6, y = env_y, label = "Inflammatory", fontface = "bold", size = 7) +
    scale_colour_manual(values = cond2_pal) +
    scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.02, 0.02)), limits = c(0, ymax)) +
    labs(title = NULL, x = NULL, y = ylab) +
    theme_prism()
  if (!is.null(br) && nrow(br)) {
    p <- p + geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = br, aes(x = (x + xend) / 2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.05, size = 5)
  }
  save_fig(p, file, w = 12, h = 6.3)
  message("OK ", file)
}

fig_dot("n_CD45",      "CD45+ total cell count",   "Fig6dot_CD45")
fig_dot("n_CD4",       "CD4+ cell count",          "Fig6dot_CD4")
fig_dot("n_CD8",       "CD8+ cell count",          "Fig6dot_CD8")
fig_dot("n_Mono_CD14", "CD14+ cell count",         "Fig6dot_Monocytes")
fig_dot("n_NK_CD16",   "CD16+ NK cell count",      "Fig6dot_NK")
fig_dot("n_B_CD19",    "CD19+ B cell count",       "Fig6dot_Bcells")
