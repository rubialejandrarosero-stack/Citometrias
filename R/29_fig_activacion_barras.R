# R/29_fig_activacion_barras.R — activation-marker figures (% of parent), filled Prism
# bars (blue = resting, orange = activated), individual donors + mean ± SD, significance
# from the linear mixed model. Marker name as title, no statistical subtitle.
source("R/fig_theme.R")
suppressPackageStartupMessages({library(lmerTest); library(emmeans); library(scales)})

cond2_lv  <- c("Resting", "Activated")
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
        axis.title.y = element_text(face = "bold", size = base + 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base + 3),
        legend.position = "none", panel.grid = element_blank())

d0 <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(environment = factor(ifelse(Condicion == "INF", "Inflammatory", "Basal"),
                              levels = c("Basal", "Inflammatory")),
         cond = factor(ifelse(Activacion == "ACT", "Activated", "Resting"), levels = cond2_lv),
         timef = factor(Tiempo), donor = factor(Donante), time = Tiempo)

fig_act <- function(pop, title, ylab, file) {
  d <- d0 %>% mutate(value = .data[[pop]]) %>% filter(!is.na(value))
  ph <- tryCatch({
    m <- lmer(value ~ environment * cond * timef + (1 | donor), data = d)
    as.data.frame(pairs(emmeans(m, ~ cond | environment * timef), adjust = "tukey"))
  }, error = function(e) NULL)

  s <- d %>% group_by(environment, time, cond) %>%
    summarise(M = mean(value), SD = sd(value), .groups = "drop") %>%
    mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])
  d <- d %>% mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])

  br <- NULL
  if (!is.null(ph)) {
    ybase <- s %>% group_by(environment, time) %>%
      summarise(yb = max(M + SD, na.rm = TRUE), .groups = "drop") %>%
      mutate(environment = as.character(environment))
    br <- ph %>% mutate(environment = as.character(environment),
                        time = as.numeric(as.character(timef)), lab = stars(p.value)) %>%
      filter(!is.na(lab)) %>% left_join(ybase, by = c("environment", "time")) %>%
      mutate(x = xbase(environment, time) + cond_off[["Resting"]],
             xend = xbase(environment, time) + cond_off[["Activated"]], y = yb * 1.06)
  }
  topbar <- max(s$M + s$SD, na.rm = TRUE)
  maxc <- if (!is.null(br) && nrow(br)) max(br$y) else topbar
  time_y <- maxc * 1.11; env_y <- time_y * 1.07; ytop <- env_y * 1.06

  axis_df <- s %>% distinct(xpos, cond) %>% mutate(lab = cond_short[as.character(cond)]) %>% arrange(xpos)
  time_df <- s %>% distinct(environment, time) %>% mutate(x = xbase(environment, time), y = time_y)

  p <- ggplot(s, aes(xpos, M, fill = cond)) +
    geom_col(width = 0.32, colour = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = M, ymax = M + SD), width = 0.14, linewidth = 0.4) +
    geom_point(data = d, aes(xpos, value), inherit.aes = FALSE,
               position = position_jitter(width = 0.06, height = 0), size = 1, colour = "black", alpha = 0.7) +
    geom_vline(xintercept = 4, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text(data = time_df, aes(x, y, label = paste0(time, " h")), inherit.aes = FALSE, fontface = "bold", size = 5) +
    annotate("text", x = 2, y = env_y, label = "Basal", fontface = "bold", size = 7) +
    annotate("text", x = 6, y = env_y, label = "Inflammatory", fontface = "bold", size = 7) +
    scale_fill_manual(values = cond2_pal) +
    scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.02)), limits = c(0, ytop)) +
    labs(title = title, x = NULL, y = ylab) +
    theme_prism()
  if (!is.null(br) && nrow(br)) {
    p <- p + geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = br, aes(x = (x + xend) / 2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.05, size = 5)
  }
  save_fig(p, file, w = 12, h = 6.3)
  message("OK ", file)
}

# NOTE: HLA-DR-based activation (CD4/CD8/CD64-HLA-DR) removed — HLA-DR does not resolve
# in this violet-heavy panel (see R/33 diagnostic: monocytes appear HLA-DR-negative while
# T cells appear HLA-DR-positive, i.e. biologically inverted). Only CD64/CD11b is reliable.
fig_act("pctPar_CD64CD11b_de_Mono", "CD64+ CD11b+", "% of CD14+ cells", "Fig6act_CD64CD11b")
