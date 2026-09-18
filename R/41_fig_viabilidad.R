# R/41_fig_viabilidad.R — viability (% live = Zombie Red-) and PD-L1 (% of spheroid CD45-)
# from the death/immune-blockade cytometry experiment. Basal | Inflammatory split,
# Control/Resting/Activated, faceted by time. Bars = mean, error = SD (2 donors), points = donors.
# Control (spheroid only) has a single value per time/environment -> no SD.
# NOTE: EpCAM is NOT plotted — its column in the source table is a broken constant (4.9%).
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(scales)})
OUT <- "output/figuras"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cond_lv <- c("Control", "Resting", "Activated"); timef_lv <- c("24", "48", "96")
GREEN  <- c(Control = "#c7e9c0", Resting = "#a1d99b", Activated = "#41ab5d")
PURPLE <- c(Control = "#dadaeb", Resting = "#bcbddc", Activated = "#807dba")

theme_ihc <- function(base = 17) theme_classic(base_size = base) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.5), axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size = base - 4),
        axis.title.y = element_text(face = "bold", size = base),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base + 3),
        strip.background = element_rect(fill = "grey92", colour = "grey70"),
        strip.text = element_text(face = "bold", size = base - 1),
        panel.border = element_rect(fill = NA, colour = "grey75", linewidth = 0.6),
        panel.spacing = unit(1.1, "lines"), legend.position = "none", panel.grid = element_blank())

xpos_fun <- function(env, cond) { b <- match(cond, cond_lv); ifelse(env == "Basal", b, b + 4) }

plot_metric <- function(csv, valcol, title, ylab, pal, out, ytop) {
  d <- read_csv(csv, show_col_types = FALSE) %>%
    rename(value = all_of(valcol)) %>%
    mutate(condition = factor(condition, cond_lv), time = factor(time, timef_lv),
           x = xpos_fun(environment, condition))
  s <- d %>% group_by(time, environment, condition, x) %>%
    summarise(M = mean(value), SD = sd(value), .groups = "drop")
  el <- merge(data.frame(env = c("Basal", "Inflammatory"), x = c(2, 6), y = ytop * 0.96),
              data.frame(time = factor(timef_lv, timef_lv)))
  p <- ggplot() +
    geom_vline(xintercept = 4, linetype = "dotted", colour = "grey80", linewidth = 0.25) +
    geom_col(data = s, aes(x, M, fill = condition), width = 0.72, colour = "grey35", linewidth = 0.3) +
    geom_errorbar(data = filter(s, !is.na(SD)), aes(x, ymin = M - SD, ymax = M + SD), width = 0.2, linewidth = 0.4) +
    geom_point(data = d, aes(x, value), size = 1.4, colour = "grey25", position = position_jitter(width = 0.12, height = 0)) +
    geom_text(data = el, aes(x, y, label = env), fontface = "bold", size = 4.4) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7), labels = rep(c("Control", "Resting", "Activated"), 2)) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.02)), limits = c(0, ytop)) +
    facet_wrap(~ time, nrow = 1, labeller = labeller(.default = function(x) paste0(x, " h"))) +
    labs(title = title, x = NULL, y = ylab) + theme_ihc()
  ggsave(file.path(OUT, out), p, width = 12, height = 5, dpi = 300)
  ggsave(file.path(OUT, sub("png$", "pdf", out)), p, width = 12, height = 5)
  message("OK ", out)
}

plot_metric("output/tidy/viabilidad.csv", "viability", "Cell viability", "% viable cells",
            GREEN, "Fig_viability.png", 108)
plot_metric("output/tidy/pdl1_flow.csv", "pdl1", "PD-L1 expression (flow)", "% PD-L1+ (of spheroid CD45-)",
            PURPLE, "Fig_pdl1_flow.png", 100)
