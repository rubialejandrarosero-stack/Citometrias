# R/25_fig_citometria_barras_composicion.R — count bar panels styled for the montage
# Fig6_composition_counts: hollow bars with thick coloured outlines (blue = resting,
# orange = activated), no statistical subtitle, population name as the only title.
# Significance from the same linear mixed model (env x activation x time + (1|donor)).
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
        axis.title = element_text(face = "bold", size = base + 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base + 3),
        legend.position = "none", panel.grid = element_blank())

d0 <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(environment = factor(ifelse(Condicion == "INF", "Inflammatory", "Basal"),
                              levels = c("Basal", "Inflammatory")),
         cond = factor(ifelse(Activacion == "ACT", "Activated", "Resting"), levels = cond2_lv),
         timef = factor(Tiempo), donor = factor(Donante), time = Tiempo)

fig_cyto_bar <- function(pop, title, file) {
  d <- d0 %>% mutate(value = .data[[pop]]) %>% filter(!is.na(value))
  res <- tryCatch({
    m <- lmer(value ~ environment * cond * timef + (1 | donor), data = d)
    list(ph = as.data.frame(pairs(emmeans(m, ~ cond | environment * timef), adjust = "tukey")))
  }, error = function(e) NULL)

  s <- d %>% group_by(environment, time, cond) %>%
    summarise(M = mean(value), SD = sd(value), .groups = "drop") %>%
    mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])
  d <- d %>% mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])

  br <- NULL
  if (!is.null(res)) {
    ybase <- s %>% group_by(environment, time) %>%
      summarise(yb = max(M + SD, na.rm = TRUE), .groups = "drop") %>%
      mutate(environment = as.character(environment))
    br <- res$ph %>%
      mutate(environment = as.character(environment), time = as.numeric(as.character(timef)),
             lab = stars(p.value)) %>% filter(!is.na(lab)) %>%
      left_join(ybase, by = c("environment", "time")) %>%
      mutate(x = xbase(environment, time) + cond_off[["Resting"]],
             xend = xbase(environment, time) + cond_off[["Activated"]], y = yb * 1.06)
  }
  topbar <- max(s$M + s$SD, na.rm = TRUE)
  maxc <- if (!is.null(br) && nrow(br)) max(br$y) else topbar
  time_y <- maxc * 1.11; env_y <- time_y * 1.07; ytop <- env_y * 1.06

  axis_df <- s %>% distinct(xpos, cond) %>% mutate(lab = cond_short[as.character(cond)]) %>% arrange(xpos)
  time_df <- s %>% distinct(environment, time) %>% mutate(x = xbase(environment, time), y = time_y)

  p <- ggplot(s, aes(xpos, M, colour = cond)) +
    geom_col(width = 0.32, fill = NA, linewidth = 1.2) +
    geom_errorbar(aes(ymin = M, ymax = M + SD), width = 0.14, linewidth = 0.5, colour = "black") +
    geom_point(data = d, aes(xpos, value), inherit.aes = FALSE,
               position = position_jitter(width = 0.06, height = 0), size = 1, colour = "black", alpha = 0.7) +
    geom_vline(xintercept = 4, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text(data = time_df, aes(x, y, label = paste0(time, " h")), inherit.aes = FALSE, fontface = "bold", size = 5) +
    annotate("text", x = 2, y = env_y, label = "Basal", fontface = "bold", size = 7) +
    annotate("text", x = 6, y = env_y, label = "Inflammatory", fontface = "bold", size = 7) +
    scale_colour_manual(values = cond2_pal) +
    scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.02)), limits = c(0, ytop)) +
    labs(title = title, x = NULL, y = "Cell count") +
    theme_prism()
  if (!is.null(br) && nrow(br)) {
    p <- p + geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4, colour = "black") +
      geom_text(data = br, aes(x = (x + xend)/2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.05, size = 5, colour = "black")
  }
  save_fig(p, file, w = 12, h = 6)
  message("OK ", file)
}

panels <- list(
  c("n_CD45",      "Total leukocytes (CD45+)", "Fig6comp_CD45"),
  c("n_CD3",       "CD3+ T cells",             "Fig6comp_CD3"),
  c("n_CD4",       "CD4+ T cells",             "Fig6comp_CD4"),
  c("n_Mono_CD14", "Macrophages (CD14+)",      "Fig6comp_Monocytes"),
  c("n_CD8",       "CD8+ T cells",             "Fig6comp_CD8"),
  c("n_NK_CD16",   "NK cells (CD16+)",         "Fig6comp_NK"))
for (x in panels) fig_cyto_bar(x[1], x[2], x[3])
