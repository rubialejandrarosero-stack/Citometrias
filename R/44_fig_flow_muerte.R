# R/44_fig_flow_muerte.R — death / immune-blockade cytometry (rectified table).
#   (1) Viable-cell counts over time: CD45- (spheroid) vs CD45+ (PBMC) as thin lines.
#   (2) PD-L1+, EpCAM+ (of CD45- spheroid) and CD3+/PD-1+ (of CD45+ T cells) as % dot plots.
# Basal | Inflammatory split, Control/Resting/Activated, points = donors (D1/D2).
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr); library(scales); library(tidyr)})
OUT <- "output/figuras"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cond_lv <- c("Control", "Resting", "Activated"); time_lv <- c("24", "48", "96")

# ---------- (1) Viable spheroid (CD45-) as % of viable cells — bars, 3 conditions ----------
COND_COL <- c("Control" = "black", "Resting" = "#4292c6", "Activated" = "#e6550d")
FILL_COL <- c("Control" = "#d9d9d9", "Resting" = "#c6dbef", "Activated" = "#fdd0a2")  # lighter fills
xoff <- c("Control" = -0.26, "Resting" = 0, "Activated" = 0.26)
cc <- read_csv("output/tidy/cd45_counts.csv", show_col_types = FALSE) %>%
  filter(time %in% c(24, 48, 96)) %>%
  tidyr::pivot_wider(names_from = population, values_from = count) %>%
  mutate(pct = 100 * `CD45- (spheroid)` / (`CD45- (spheroid)` + `CD45+ (PBMC)`),
         condition = factor(condition, cond_lv), environment = factor(environment, c("Basal", "Inflammatory")),
         xd = match(as.character(time), time_lv) + xoff[as.character(condition)])
cs <- cc %>% group_by(time, environment, condition, xd) %>%
  summarise(M = mean(pct, na.rm = TRUE), SD = sd(pct, na.rm = TRUE), .groups = "drop")

pline <- ggplot() +
  geom_col(data = cs, aes(xd, M, fill = condition, colour = condition), width = 0.24, linewidth = 0.6) +
  geom_errorbar(data = filter(cs, !is.na(SD)), aes(xd, ymin = M, ymax = M + SD, colour = condition),
                width = 0.10, linewidth = 0.4) +
  geom_point(data = cc, aes(xd, pct, colour = condition), size = 1.0, alpha = 0.5,
             position = position_jitter(width = 0.04, height = 0)) +
  scale_colour_manual(values = COND_COL) + scale_fill_manual(values = FILL_COL) +
  scale_x_continuous(breaks = 1:3, labels = c("24", "48", "96")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, NA)) +
  facet_wrap(~ environment, nrow = 1) +
  labs(x = "Time (h)", y = "Viable CD45⁻ cells (%)", colour = NULL, fill = NULL) +
  theme_classic(base_size = 15) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.5), axis.text = element_text(colour = "black"),
        axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        strip.background = element_rect(fill = "grey92", colour = "grey70"), strip.text = element_text(face = "bold"),
        panel.border = element_rect(fill = NA, colour = "grey75", linewidth = 0.6), panel.spacing = unit(0.9, "lines"),
        legend.position = "top", panel.grid = element_blank())
# pline is saved later, panel-aligned with the EpCAM graph (see end of file)

# ---------- (2) % marker dot plots — morphology style: single graph, Basal | Inflammatory split,
#             each condition grouped with its 3 timepoints (colour = time), mean +/- SD ----------
time_off <- c("24" = -0.26, "48" = 0, "96" = 0.26)
TIME_COL <- c("24" = "#bdd7e7", "48" = "#4292c6", "96" = "#08306b")

dotplot <- function(csv, valcol, title, ylab, out) {
  d <- read_csv(csv, show_col_types = FALSE) %>%
    filter(time %in% c(24, 48, 96)) %>%
    rename(value = all_of(valcol)) %>%
    mutate(condition = factor(condition, cond_lv), timef = factor(time, time_lv),
           x = match(condition, cond_lv) + ifelse(environment == "Basal", 0, 4) +
               time_off[as.character(time)])
  s <- d %>% group_by(environment, condition, timef, x) %>%
    summarise(M = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE), .groups = "drop")
  ytop <- max(d$value, na.rm = TRUE)
  el <- data.frame(x = c(2, 6), y = ytop * 1.1, lab = c("Basal", "Inflammatory"))
  p <- ggplot() +
    geom_vline(xintercept = 4, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_point(data = d, aes(x, value, colour = timef), shape = 1, size = 1.9, stroke = 0.9,
               alpha = 0.85, position = position_jitter(width = 0.055, height = 0)) +
    geom_errorbar(data = filter(s, !is.na(SD)), aes(x, ymin = M - SD, ymax = M + SD, colour = timef),
                  width = 0.13, linewidth = 0.55) +
    geom_errorbar(data = s, aes(x, ymin = M, ymax = M, colour = timef), width = 0.22, linewidth = 0.9) +
    geom_text(data = el, aes(x, y, label = lab), fontface = "bold", size = 5) +
    scale_colour_manual(values = TIME_COL, labels = paste0(time_lv, " h"), name = NULL) +
    scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7), labels = rep(cond_lv, 2)) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0.02, 0.15))) +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(base_size = 15) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5), axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.title.y = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          legend.position = "top", panel.grid = element_blank())
  ggsave(file.path(OUT, out), p, width = 8.5, height = 5, dpi = 300)
  ggsave(file.path(OUT, sub("png$", "pdf", out)), p, width = 8.5, height = 5)
  message("OK ", out)
}

# hollow-bar variant (same layout/colours as dotplot but unfilled bars = mean, whisker = +SD)
barplot_hollow <- function(csv, valcol, title, ylab, out) {
  d <- read_csv(csv, show_col_types = FALSE) %>%
    filter(time %in% c(24, 48, 96)) %>%
    rename(value = all_of(valcol)) %>%
    mutate(condition = factor(condition, cond_lv), timef = factor(time, time_lv),
           x = match(condition, cond_lv) + ifelse(environment == "Basal", 0, 4) +
               time_off[as.character(time)])
  s <- d %>% group_by(environment, condition, timef, x) %>%
    summarise(M = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE), .groups = "drop")
  ytop <- max(d$value, na.rm = TRUE)
  el <- data.frame(x = c(2, 6), y = ytop * 1.1, lab = c("Basal", "Inflammatory"))
  p <- ggplot() +
    geom_vline(xintercept = 4, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_col(data = s, aes(x, M, colour = timef), fill = NA, linewidth = 0.9, width = 0.24) +
    geom_errorbar(data = filter(s, !is.na(SD)), aes(x, ymin = M, ymax = M + SD, colour = timef),
                  width = 0.10, linewidth = 0.55) +
    geom_point(data = d, aes(x, value, colour = timef), shape = 1, size = 1.5, stroke = 0.8,
               alpha = 0.8, position = position_jitter(width = 0.04, height = 0)) +
    geom_text(data = el, aes(x, y, label = lab), fontface = "bold", size = 5) +
    scale_colour_manual(values = TIME_COL, labels = paste0(time_lv, " h"), name = NULL) +
    scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7), labels = rep(cond_lv, 2)) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.15))) +
    labs(title = title, x = NULL, y = ylab) +
    theme_classic(base_size = 15) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5), axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.title.y = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          legend.position = "top", panel.grid = element_blank())
  ggsave(file.path(OUT, out), p, width = 8.5, height = 5, dpi = 300)
  ggsave(file.path(OUT, sub("png$", "pdf", out)), p, width = 8.5, height = 5)
  message("OK ", out)
}

# line variant: 3 condition lines joining the timepoints, Basal | Inflammatory facets (as CD45 lines)
lineplot <- function(csv, valcol, title, ylab, out, show_legend = TRUE, save = TRUE) {
  d <- read_csv(csv, show_col_types = FALSE) %>%
    filter(time %in% c(24, 48, 96)) %>%
    rename(value = all_of(valcol)) %>%
    mutate(condition = factor(condition, cond_lv), environment = factor(environment, c("Basal", "Inflammatory")))
  s <- d %>% group_by(time, environment, condition) %>%
    summarise(M = mean(value, na.rm = TRUE), SD = sd(value, na.rm = TRUE), .groups = "drop")
  p <- ggplot(s, aes(time, M, colour = condition)) +
    geom_errorbar(aes(ymin = pmax(M - SD, 0), ymax = M + SD), width = 3, linewidth = 0.35, na.rm = TRUE) +
    geom_line(linewidth = 0.7) +
    geom_point(data = d, aes(time, value), size = 1.1, alpha = 0.5,
               position = position_jitter(width = 1.2, height = 0)) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = COND_COL) +
    scale_x_continuous(breaks = c(24, 48, 96)) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.08)), limits = c(0, NA)) +
    facet_wrap(~ environment, nrow = 1) +
    labs(x = "Time (h)", y = ylab, colour = NULL, title = title) +
    theme_classic(base_size = 15) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5), axis.text = element_text(colour = "black"),
          axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          strip.background = element_rect(fill = "grey92", colour = "grey70"), strip.text = element_text(face = "bold"),
          panel.border = element_rect(fill = NA, colour = "grey75", linewidth = 0.6), panel.spacing = unit(0.9, "lines"),
          legend.position = "top", panel.grid = element_blank())
  if (!show_legend)   # keep the legend row (so gtable rows match) but render it invisible (white)
    p <- p + guides(colour = guide_legend(override.aes = list(colour = "white"))) +
      theme(legend.text = element_text(colour = "white"), legend.key = element_rect(fill = "white", colour = "white"))
  if (save) {
    ggsave(file.path(OUT, out), p, width = 8, height = 5.2, dpi = 300)
    ggsave(file.path(OUT, sub("png$", "pdf", out)), p, width = 8, height = 5.2)
    message("OK ", out)
  }
  invisible(p)
}

save_grob <- function(g, file) {
  png(file.path(OUT, file), width = 8, height = 5.2, units = "in", res = 300, bg = "white")
  grid::grid.newpage(); grid::grid.draw(g); dev.off()
}
align_pair <- function(p1, p2) {   # force identical panel size (widths + heights) via gtable
  g1 <- ggplotGrob(p1); g2 <- ggplotGrob(p2)
  mw <- grid::unit.pmax(g1$widths, g2$widths);   g1$widths  <- mw; g2$widths  <- mw
  mh <- grid::unit.pmax(g1$heights, g2$heights); g1$heights <- mh; g2$heights <- mh
  list(g1, g2)
}

# Figure 4 graphs (CD45 viability + EpCAM): no titles, legend only on CD45, identical panel size
epcam_p <- lineplot("output/tidy/epcam_flow.csv", "epcam", NULL, "% EpCAM-positive cells",
                    "Fig_epcam_flow.png", show_legend = FALSE, save = FALSE)
g <- align_pair(pline, epcam_p)
save_grob(g[[1]], "Fig_cd45_viable_lines.png"); save_grob(g[[2]], "Fig_epcam_flow.png")
message("OK CD45 + EpCAM (panel-aligned)")

# Figure 9 graphs (PD-L1 + PD-1 flow): no titles, legend only on PD-L1, identical panel size
pdl1_p <- lineplot("output/tidy/pdl1_flow.csv", "pdl1", NULL, "% PD-L1-positive cells",
                   "Fig_pdl1_flow.png", show_legend = TRUE, save = FALSE)
pd1_p  <- lineplot("output/tidy/pd1_flow.csv", "pd1", NULL, "% PD-1-positive cells",
                   "Fig_pd1_flow.png", show_legend = FALSE, save = FALSE)
g <- align_pair(pdl1_p, pd1_p)
save_grob(g[[1]], "Fig_pdl1_flow.png"); save_grob(g[[2]], "Fig_pd1_flow.png")
message("OK PD-L1 + PD-1 (panel-aligned)")
