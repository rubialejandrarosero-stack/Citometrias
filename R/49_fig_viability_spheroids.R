# R/49_fig_viability_spheroids.R — spheroid viability (Live/Dead of Singlets) from the FlowJo
# gating export "Viabilidad esferoides 08 sept.xls". Percent = of Singlets (parent gate) and
# absolute counts (#Cells). 2 donors (D1, D2), 3 timepoints (24/48/96 h), Basal | Inflammatory x
# Resting | Activated, plus Control (ESF SOLO, no PBMC). Same layout/palette as R/48.
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(grid)
})
OUT <- "output/figuras"
XLS <- "/mnt/c/Users/57319/OneDrive/Escritorio/Viabilidad esferoides 08 sept.xls"

cond_lv  <- c("Control", "Resting", "Activated")
cond_pal <- c(Control = "black", Resting = "#4292c6", Activated = "#e31a1c")
cond_lab <- c(Control = "Control", Resting = "Resting PBMC", Activated = "Activated PBMC")
cond_off <- c(Control = -0.32, Resting = 0, Activated = 0.32)
env_lab  <- c(Basal = "Basal", Inflammatory = "TNF α, IL-α, IL-1 β")

# ---- parse the FlowJo export ----
raw <- read_excel(XLS, sheet = "Sheet0", col_names = TRUE)[, 1:4]
names(raw) <- c("depth", "name", "stat", "ncells")
raw <- raw %>% filter(!is.na(name)) %>%
  mutate(stat = as.numeric(stat), ncells = as.numeric(ncells))

# drop the two non-experimental pseudo-samples (SIN MARCAR = unstained control,
# Tumorales.fcs = FlowJo's concatenated pool, both technical, not experimental conditions)
raw <- raw %>% filter(!str_detect(name, "^SIN MARCAR|^Tumorales\\.fcs"))

parse_sample <- function(fname) {
  base <- str_remove(fname, "\\.fcs$")
  time <- as.numeric(str_extract(base, "(?<=_)\\d+$"))
  head <- str_remove(base, "_\\d+$")
  if (str_detect(head, "^ESF SOLO")) {
    environment <- ifelse(str_detect(head, "NO INF"), "Basal", "Inflammatory")
    condition <- "Control"; donor <- NA_character_
  } else {
    environment <- ifelse(str_detect(head, "^NO INF"), "Basal", "Inflammatory")
    condition <- ifelse(str_detect(head, "NO ACT"), "Resting", "Activated")
    donor <- str_extract(head, "D\\d+")
  }
  tibble(environment = environment, condition = condition, donor = donor, time = time)
}

raw <- raw %>% mutate(sample = str_extract(name, "^[^/]+"), path = str_remove(name, "^[^/]+/?"))
meta <- raw %>% distinct(sample) %>% rowwise() %>% mutate(parse_sample(sample)) %>% ungroup()

d <- raw %>% left_join(meta, by = "sample") %>%
  mutate(condition = factor(condition, cond_lv), environment = factor(environment, c("Basal", "Inflammatory")))

nodes <- c(live = "cells/Singlets/Live", dead = "cells/Singlets/Dead")

get_node <- function(tag) {
  d %>% filter(path == nodes[[tag]]) %>%
    transmute(environment, condition, donor, time, pct = stat, n = ncells, pop = tag)
}
pops <- bind_rows(lapply(names(nodes), get_node))

readr::write_csv(pops, file.path("output/tidy", "viability_populations.csv"))

summ <- pops %>% group_by(pop, environment, condition, time) %>%
  summarise(M_pct = mean(pct, na.rm = TRUE), SD_pct = sd(pct, na.rm = TRUE),
            M_n = mean(n, na.rm = TRUE), SD_n = sd(n, na.rm = TRUE), .groups = "drop")

tidx <- function(t) match(t, c(24, 48, 96)) * 1.4
BLOCK_OFF <- 2.6
INF_SHIFT <- tidx(96) + BLOCK_OFF - tidx(24)
xbase <- function(env, t) tidx(t) + ifelse(env == "Basal", 0, INF_SHIFT)

panel <- function(df_pts, df_sum, y_pts, y_m, y_sd, ylab, title, bar_width = 0.28, offsets = cond_off) {
  dd <- df_pts %>% mutate(x = xbase(environment, time) + offsets[as.character(condition)],
                          val = .data[[y_pts]])
  ss <- df_sum %>% mutate(x = xbase(environment, time) + offsets[as.character(condition)],
                          M = .data[[y_m]], SD = .data[[y_sd]])
  ggplot() +
    geom_col(data = ss, aes(x, M, colour = condition), fill = NA, width = bar_width, linewidth = 1.0) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + replace_na(SD, 0), colour = condition),
                  width = bar_width * 0.46, linewidth = 0.5) +
    geom_point(data = dd, aes(x, val, colour = condition), size = 1.6, alpha = 0.85) +
    geom_vline(xintercept = (xbase("Basal", 96) + xbase("Inflammatory", 24)) / 2,
              colour = "black", linewidth = 1.0) +
    scale_colour_manual(values = cond_pal, labels = cond_lab, name = NULL) +
    scale_x_continuous(breaks = c(xbase("Basal", c(24, 48, 96)), xbase("Inflammatory", c(24, 48, 96))),
                       labels = rep(c("24", "48", "96"), 2),
                       limits = c(tidx(24) - 0.6, xbase("Inflammatory", 96) + 0.6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Time (h)", y = ylab, title = title) +
    theme_classic(base_size = 15) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.x = element_blank(),
          axis.text = element_text(colour = "black", size = 13),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.title.y = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
          legend.position = "top", legend.direction = "horizontal",
          legend.text = element_text(size = 13), panel.grid = element_blank())
}

pop_titles <- c(live = "Spheroid viability — live (Singlets)", dead = "Spheroid viability — dead (Singlets)")

render_grid <- function(measure, ylab_pct, ylab_n, outfile,
                        pct_breaks = NULL, pct_max = NULL, n_breaks = NULL, n_max = NULL,
                        bar_width = 0.28, drop_control = FALSE) {
  pts <- pops %>% filter(pop == measure)
  sm  <- summ  %>% filter(pop == measure)
  offs <- cond_off
  if (drop_control) {
    pts <- filter(pts, condition != "Control") %>% mutate(condition = droplevels(condition))
    sm  <- filter(sm,  condition != "Control") %>% mutate(condition = droplevels(condition))
    offs <- c(Resting = -0.19, Activated = 0.19)
  }
  p_pct <- panel(pts, sm, "pct", "M_pct", "SD_pct", ylab_pct, pop_titles[[measure]], bar_width, offs)
  p_n   <- panel(pts, sm, "n",   "M_n",   "SD_n",   ylab_n,   pop_titles[[measure]], bar_width, offs)

  p_pct <- p_pct + if (is.null(pct_breaks)) {
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)))
  } else {
    scale_y_continuous(breaks = pct_breaks, limits = c(0, pct_max), expand = expansion(mult = c(0, 0.02)))
  }
  p_n <- p_n + if (is.null(n_breaks)) {
    scale_y_continuous(breaks = scales::breaks_pretty(n = 6), expand = expansion(mult = c(0, 0.08)))
  } else {
    scale_y_continuous(breaks = n_breaks, limits = c(0, n_max), expand = expansion(mult = c(0, 0.02)))
  }

  gl <- ggplotGrob(p_pct + theme(legend.position = "bottom", legend.direction = "horizontal"))
  legend <- gl$grobs[[which(vapply(gl$grobs, function(x) x$name, "") == "guide-box")]]

  png(file.path(OUT, outfile), width = 11.5, height = 5.8, units = "in", res = 300, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2, heights = unit.c(unit(1, "null"), unit(2.2, "lines")))))
  prt <- function(p, r, c) { pushViewport(viewport(layout.pos.row = r, layout.pos.col = c)); print(p, newpage = FALSE); popViewport() }
  prt(p_pct + theme(legend.position = "none"), 1, 1)
  prt(p_n   + theme(legend.position = "none"), 1, 2)
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2)); grid.draw(legend); popViewport()
  popViewport(); dev.off()
  message("OK ", outfile)
}

render_grid("live", "% of Singlets", "Absolute count", "Fig_spheroid_viability_live.png",
           pct_breaks = seq(0, 60, 10), pct_max = 60, n_breaks = seq(0, 2500, 500), n_max = 2500)
render_grid("dead", "% of Singlets", "Absolute count", "Fig_spheroid_viability_dead.png")
