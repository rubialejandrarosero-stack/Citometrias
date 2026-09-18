# R/43_fig_cytokines.R — key cytokines from the Luminex secretome (pg/mL): CD25, CXCL10, CXCL9,
# IFN-g, IL-10, IL-17, FGF-2, GM-CSF. Each cytokine is an INDEPENDENT panel with its own plain
# title and legend, broken down by its 3 timepoints (24/48/96 h), Basal | Inflammatory groups
# (thick divider + in-panel "No inflamatorio"/"Inflamatorio" label), Control/Resting/Activated as
# dodged HOLLOW bars (colour = condition, black/blue/red) + donor points + SD. Arranged 2 x 4.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales); library(grid); library(tidyr)
})
OUT <- "output/figuras"
cond_lv  <- c("Control", "Resting", "Activated")
cond_pal <- c(Control = "black", Resting = "#4292c6", Activated = "#e31a1c")
cond_lab <- c("Control", "Resting PBMC", "Activated PBMC")
cond_off <- c(Control = -0.24, Resting = 0, Activated = 0.24)
env_lab  <- c(Basal = "Basal", Inflammatory = "TNF α, IL-α, IL-1 β")
tidx <- function(t) match(t, c(24, 48, 96))
xbase <- function(env, t) ifelse(env == "Basal", tidx(t), tidx(t) + 4)
KEY <- c("CD25" = "CD25", "CXCL10/IP10/CRG2" = "CXCL10", "CXCL-9" = "CXCL9",
         "IFN-gamma" = "IFN-γ", "IL-10" = "IL-10", "IL-17" = "IL-17",
         "FGF basic/ FGF2/bFGF" = "FGF-2", "GM-CSF" = "GM-CSF")

d0 <- read_csv("output/tidy/luminex_tidy.csv", show_col_types = FALSE) %>%
  filter(cytokine %in% names(KEY)) %>%
  mutate(condition = factor(ifelse(pbmc %in% c("None", NA), "Control", pbmc), cond_lv),
         cyt = factor(KEY[cytokine], levels = unname(KEY)),
         environment = factor(environment, c("Basal", "Inflammatory")), timef = factor(time)) %>%
  group_by(cytokine, environment, condition, time) %>% mutate(rep_id = row_number()) %>% ungroup() %>%
  mutate(donor2 = factor(ifelse(condition == "Control", paste0("C", rep_id), as.character(donor))),
         x = xbase(environment, time) + cond_off[as.character(condition)])
s0 <- d0 %>% group_by(cyt, environment, condition, time, x) %>%
  summarise(M = mean(value), SD = sd(value), .groups = "drop")

panel <- function(cy) {
  dd <- filter(d0, cyt == cy); ss <- filter(s0, cyt == cy)
  bar_top <- max(ss$M + replace_na(ss$SD, 0), dd$value, na.rm = TRUE)
  ytop <- bar_top * 1.18

  # in-panel environment labels ("No inflamatorio" / "Inflamatorio"), centred over each half, near the top
  env_labs <- tibble(x = c(2, 6), lab = c("No inflamatorio", "Inflamatorio"), y = ytop * 0.98)

  p <- ggplot() +
    geom_col(data = ss, aes(x, M, colour = condition), fill = NA, width = 0.20, linewidth = 1.8) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + SD, colour = condition),
                  width = 0.12, linewidth = 1.1) +
    geom_point(data = dd, aes(x, value, colour = condition), size = 1.1, alpha = 0.8,
               position = position_jitter(width = 0.05, height = 0)) +
    geom_vline(xintercept = 4, colour = "black", linewidth = 1.0) +      # thick divider Basal | cocktail
    geom_text(data = env_labs, aes(x, y, label = lab), inherit.aes = FALSE,
              size = 12, fontface = "italic", colour = "grey30", vjust = 1) +
    scale_colour_manual(values = cond_pal, labels = cond_lab, name = NULL) +
    scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7), labels = rep(c("24 h", "48 h", "96 h"), 2),
                       limits = c(0.0, 7.85)) +
    scale_y_continuous(labels = label_comma(), limits = c(0, ytop),
                       expand = expansion(mult = c(0, 0.02))) +      # hard floor at 0 -> bars sit on the axis
    labs(x = NULL, y = "Concentration (pg/mL)", title = as.character(cy)) +
    theme_classic(base_size = 40) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.x = element_blank(),
          axis.text = element_text(colour = "black", size = 40),
          axis.text.x = element_text(size = 34, angle = 45, hjust = 1),
          axis.title.y = element_text(face = "bold", size = 40),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 40),
          legend.position = "top", legend.direction = "horizontal",
          legend.text = element_text(size = 40), legend.key.size = unit(1.8, "lines"),
          legend.margin = margin(1, 1, 1, 1), panel.grid = element_blank())

  p
}

cys <- unname(KEY)
png(file.path(OUT, "Fig_cytokines_key.png"), width = 34.0, height = 42.0, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(7, 2,
             heights = unit.c(unit(1, "null"), unit(2.4, "lines"), unit(1, "null"), unit(2.4, "lines"),
                              unit(1, "null"), unit(2.4, "lines"), unit(1, "null")))))
prt <- function(p, r, c) {
  pushViewport(viewport(layout.pos.row = r, layout.pos.col = c, clip = "off"))  # let below-axis daggers bleed out
  print(p, newpage = FALSE); popViewport()
}
panel_rows <- c(1, 3, 5, 7)
for (i in seq_along(cys)) {
  r <- panel_rows[(i - 1) %/% 2 + 1]
  prt(panel(cys[i]), r, (i - 1) %% 2 + 1)
}
popViewport(); dev.off()
message("OK Fig_cytokines_key.png")

# Individual per-cytokine files (same panel design/size as one cell of the combined grid)
IND_DIR <- file.path(OUT, "cytokines_individual")
dir.create(IND_DIR, showWarnings = FALSE, recursive = TRUE)
slug <- function(s) {
  s <- gsub("α", "alpha", s); s <- gsub("β", "beta", s); s <- gsub("γ", "gamma", s)
  gsub("[^A-Za-z0-9]+", "", s)
}
for (i in seq_along(cys)) {
  fn <- file.path(IND_DIR, paste0("Fig_cytokine_", slug(cys[i]), ".png"))
  ggsave(fn, panel(cys[i]), width = 17.0, height = 12.0, units = "in", dpi = 300, bg = "white")
  message("OK ", fn)
}
