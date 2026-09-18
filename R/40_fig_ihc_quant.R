# R/40_fig_ihc_quant.R — IHC quantification bar charts (pastel blues), laid out like the
# IHC image compositions (Basal | Inflammatory split). Bars = mean of biological replicates,
# error bars = SD, individual replicate points overlaid. No inner gridlines / no dashed dividers
# (facet panel borders separate the timepoints). CD14/CD3 uses a split (two-section) y-axis so
# the small bars remain readable next to the large ones.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales); library(grid)
})
OUT <- "output/figuras"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
BLUE1 <- "#c6dbef"; BLUE2 <- "#6baed6"; SINGLE <- "#9ecae1"     # pastel blues (CD14/CD3)
cond_pal <- c(Control = "#cfe3f2", Resting = "#93c4de", Activated = "#4a97c9")  # pastel blue per condition
cond_lv <- c("Control", "Resting", "Activated"); timef_lv <- c("24", "48", "96")

theme_ihc <- function(base = 17) theme_classic(base_size = base) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base - 4),
        axis.title.y = element_text(face = "bold", size = base),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base + 3),
        strip.background = element_rect(fill = "grey92", colour = "grey70"),
        strip.text = element_text(face = "bold", size = base - 1),
        panel.border = element_rect(fill = NA, colour = "grey75", linewidth = 0.6),
        panel.spacing = unit(1.1, "lines"),
        legend.position = "top", legend.title = element_blank(),
        panel.grid = element_blank())

xpos_fun <- function(env, cond, conds) {
  base <- match(cond, conds); nc <- length(conds)
  ifelse(env == "Basal", base, base + nc + 1)
}
env_labels <- function(nc) data.frame(env = c("Basal", "Inflammatory"),
                                       x = c((1 + nc) / 2, (nc + 2 + 2 * nc + 1) / 2))

# ---- single-panel charts (Ki-67 48 h, PD-L1 48 h): mean + SD + replicate points ----
single_bar <- function(d, nc, title, file, w, h) {
  conds <- cond_lv[seq_len(nc)]
  d <- d %>% mutate(condition = factor(condition, conds), x = xpos_fun(environment, condition, conds))
  s <- d %>% group_by(environment, condition, x) %>% summarise(M = mean(value), SD = sd(value), .groups = "drop")
  ytop <- max(s$M + s$SD, na.rm = TRUE) * 1.20
  el <- env_labels(nc) %>% mutate(y = max(s$M + s$SD, na.rm = TRUE) * 1.12)
  p <- ggplot() +
    geom_vline(xintercept = nc + 1, linetype = "dotted", colour = "grey80", linewidth = 0.25) +
    geom_col(data = s, aes(x, M, fill = condition), width = 0.72, colour = "grey35", linewidth = 0.3) +
    geom_errorbar(data = s, aes(x, ymin = M, ymax = M + SD), width = 0.2, linewidth = 0.4) +
    geom_point(data = d, aes(x, value), size = 1.2, colour = "grey25",
               position = position_jitter(width = 0.12, height = 0)) +
    geom_text(data = el, aes(x, y, label = env), fontface = "bold", size = 4.6) +
    scale_fill_manual(values = cond_pal) +
    scale_x_continuous(breaks = s$x, labels = as.character(s$condition)) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.02)), limits = c(0, ytop)) +
    labs(title = title, x = NULL, y = "% positive cells") + theme_ihc() +
    theme(legend.position = "none")
  ggsave(file.path(OUT, file), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, sub("png$", "pdf", file)), p, width = w, height = h)
  message("OK ", file)
}

single_bar(read_csv("output/tidy/ihc_ki67.csv", show_col_types = FALSE) %>% filter(time == 48),
           3, "Ki-67 proliferation (48 h)", "Fig_quant_Ki67.png", 6.5, 5.2)
single_bar(read_csv("output/tidy/ihc_pdl1.csv", show_col_types = FALSE),
           3, "PD-L1 expression (48 h)", "Fig_quant_PDL1.png", 6.5, 5.2)

# ---- CD14 / CD3: one hollow-bar chart per environment, TIME on the x-axis (shows the time course) ----
cd_conds <- c("Resting", "Activated")
env_bars <- function(envname, out) {
  d <- read_csv("output/tidy/ihc_cd14cd3.csv", show_col_types = FALSE) %>%
    filter(environment == envname) %>%
    mutate(condition = factor(condition, cd_conds), marker = factor(marker, c("CD14", "CD3")),
           ti = match(as.character(time), timef_lv),
           xd = ti + ifelse(marker == "CD14", -0.18, 0.18))
  s <- d %>% group_by(condition, marker, time, ti, xd) %>%
    summarise(M = mean(value), SD = sd(value), .groups = "drop")
  p <- ggplot() +
    geom_col(data = s, aes(xd, M, colour = marker), fill = NA, width = 0.32, linewidth = 1.4) +
    geom_errorbar(data = s, aes(xd, ymin = M, ymax = M + SD, colour = marker), width = 0.12, linewidth = 0.4) +
    geom_point(data = d, aes(xd, value), size = 0.7, colour = "grey25",
               position = position_jitter(width = 0.04, height = 0)) +
    scale_colour_manual(values = c(CD14 = "#2171b5", CD3 = "#cb181d"), name = NULL) +
    scale_x_continuous(breaks = 1:3, labels = paste0(timef_lv, " h")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, NA)) +
    facet_wrap(~ condition, nrow = 1) +
    labs(title = envname, x = "Time (h)", y = "% positive cells") +
    theme_ihc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  ggsave(file.path(OUT, out), p, width = 7, height = 4.6, dpi = 300, bg = "white")
  ggsave(file.path(OUT, sub("png$", "pdf", out)), p, width = 7, height = 4.6, bg = "white")
  message("OK ", out)
}
env_bars("Basal", "Fig_quant_CD14CD3_basal.png")
env_bars("Inflammatory", "Fig_quant_CD14CD3_inflammatory.png")
