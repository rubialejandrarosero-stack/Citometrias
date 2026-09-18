# R/50_fig_mfi_pd1_epcam.R — MFI plots for PD-1 (Comp-APC-A, of CD3+ PD-1+ live PBMC within the
# spheroid) and EpCAM (Comp-Alexa Fluor 700-A, of live/total EpCAM+ tumour cells), from two
# separate FlowJo exports. Same 2-donor / 3-timepoint / Basal|Inflammatory layout as R/48-49.
# PD-1 MFI: only Resting/Activated have a PD-1+ CD3+ PBMC population (Control has no PBMC).
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(grid)
})
OUT <- "output/figuras"
PD1_XLS   <- "/mnt/c/Users/57319/OneDrive/Escritorio/08-09-2026.wsp PD-1+.xls"
EPCAM_XLS <- "/mnt/c/Users/57319/OneDrive/Documentos/08-09-2026.wsp EPCAM MFI.xls"

cond_lv  <- c("Control", "Resting", "Activated")
cond_pal <- c(Control = "black", Resting = "#4292c6", Activated = "#e31a1c")
cond_lab <- c(Control = "Control", Resting = "Resting PBMC", Activated = "Activated PBMC")
cond_off <- c(Control = -0.32, Resting = 0, Activated = 0.32)
cond_off2 <- c(Resting = -0.19, Activated = 0.19)

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

tidx <- function(t) match(t, c(24, 48, 96)) * 1.4
BLOCK_OFF <- 2.6
INF_SHIFT <- tidx(96) + BLOCK_OFF - tidx(24)
xbase <- function(env, t) tidx(t) + ifelse(env == "Basal", 0, INF_SHIFT)

panel <- function(df_pts, df_sum, ylab, title, bar_width = 0.28, offsets = cond_off, ymax = NULL) {
  dd <- df_pts %>% mutate(x = xbase(environment, time) + offsets[as.character(condition)])
  ss <- df_sum %>% mutate(x = xbase(environment, time) + offsets[as.character(condition)])
  ggplot() +
    geom_col(data = ss, aes(x, M, colour = condition), fill = NA, width = bar_width, linewidth = 1.0) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + replace_na(SD, 0), colour = condition),
                  width = bar_width * 0.46, linewidth = 0.5) +
    geom_point(data = dd, aes(x, mfi, colour = condition), size = 1.6, alpha = 0.85) +
    geom_vline(xintercept = (xbase("Basal", 96) + xbase("Inflammatory", 24)) / 2,
              colour = "black", linewidth = 1.0) +
    scale_colour_manual(values = cond_pal, labels = cond_lab, name = NULL) +
    scale_x_continuous(breaks = c(xbase("Basal", c(24, 48, 96)), xbase("Inflammatory", c(24, 48, 96))),
                       labels = rep(c("24 h", "48 h", "96 h"), 2),
                       limits = c(tidx(24) - 0.6, xbase("Inflammatory", 96) + 0.6)) +
    { if (is.null(ymax)) scale_y_continuous(expand = expansion(mult = c(0, 0.08)))
      else scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) } +
    labs(x = NULL, y = ylab, title = title) +
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

render_single <- function(pts, sm, ylab, title, outfile, offs = cond_off, ymax = NULL) {
  p <- panel(pts, sm, ylab, title, offsets = offs, ymax = ymax)
  gl <- ggplotGrob(p + theme(legend.position = "bottom", legend.direction = "horizontal"))
  legend <- gl$grobs[[which(vapply(gl$grobs, function(x) x$name, "") == "guide-box")]]
  png(file.path(OUT, outfile), width = 7.0, height = 5.8, units = "in", res = 300, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit.c(unit(1, "null"), unit(2.2, "lines")))))
  pushViewport(viewport(layout.pos.row = 1)); print(p + theme(legend.position = "none"), newpage = FALSE); popViewport()
  pushViewport(viewport(layout.pos.row = 2)); grid.draw(legend); popViewport()
  popViewport(); dev.off()
  message("OK ", outfile)
}

# ============================ PD-1 MFI ============================
pd1_raw <- read_excel(PD1_XLS, sheet = "Sheet0", col_names = TRUE)[, 1:4]
names(pd1_raw) <- c("depth", "name", "stat", "ncells")

# the "Median" row's Name sometimes carries just "Median : Comp-APC-A = <val>" and sometimes
# "<sample prefix>.Median : Comp-APC-A = <val>" (FlowJo export format has varied between saves) ->
# always resolve the sample from the nearest preceding sample-header row (depth is NA), not from
# any prefix embedded in the Median row's own text. The MFI value itself is read from the `stat`
# column, NOT parsed from the Name text: some rows' Name text carries a stale/mismatched number
# (confirmed: two rows both showed "= 4220.45" in Name while their `stat` values were the correct,
# distinct 563 and 422) -- `stat` is the reliable per-row value.
med_idx <- which(str_detect(pd1_raw$name, "Median : Comp-APC-A") & !is.na(pd1_raw$stat) & pd1_raw$stat != "n/a")
sample_idx <- which(is.na(pd1_raw$depth) & !is.na(pd1_raw$name))

sample_for_row <- function(i) max(sample_idx[sample_idx < i])
pd1 <- tibble(
  sample = pd1_raw$name[vapply(med_idx, sample_for_row, integer(1))],
  mfi = as.numeric(pd1_raw$stat[med_idx])
) %>%
  filter(!str_detect(sample, "^SIN MARCAR|^Tumorales\\.fcs")) %>%
  rowwise() %>% mutate(parse_sample(sample)) %>% ungroup() %>%
  filter(condition != "Control") %>%                       # Control has no PBMC -> no PD-1+ CD3+ node
  mutate(condition = factor(condition, cond_lv) %>% droplevels(),
         environment = factor(environment, c("Basal", "Inflammatory")))

readr::write_csv(pd1, file.path("output/tidy", "pd1_mfi.csv"))

pd1_summ <- pd1 %>% group_by(environment, condition, time) %>%
  summarise(M = mean(mfi, na.rm = TRUE), SD = sd(mfi, na.rm = TRUE), .groups = "drop")

render_single(pd1, pd1_summ, "PD-1 MFI (Comp-APC-A)", "PD-1 MFI (CD3+ PD-1+ live PBMC)",
             "Fig_pd1_mfi_new.png", offs = cond_off2)

# ============================ EpCAM MFI ============================
epcam_raw <- read_excel(EPCAM_XLS, sheet = "Sheet0", col_names = TRUE)
names(epcam_raw) <- c("sample", "mfi")
epcam <- epcam_raw %>%
  filter(!sample %in% c("SIN MARCAR.fcs", "Tumorales.fcs", "Mean", "SD"), !is.na(mfi)) %>%
  rowwise() %>% mutate(parse_sample(sample)) %>% ungroup() %>%
  mutate(condition = factor(condition, cond_lv), environment = factor(environment, c("Basal", "Inflammatory")))

readr::write_csv(epcam, file.path("output/tidy", "epcam_mfi_new.csv"))

epcam_summ <- epcam %>% group_by(environment, condition, time) %>%
  summarise(M = mean(mfi, na.rm = TRUE), SD = sd(mfi, na.rm = TRUE), .groups = "drop")

render_single(epcam, epcam_summ, "EpCAM MFI (Comp-Alexa Fluor 700-A)", "EpCAM MFI (Tumour cells)",
             "Fig_epcam_mfi_new.png", offs = cond_off, ymax = 4000)

# ============================ PD-1 % and absolute count (of CD3+ live PBMC) ============================
# replaces the pd1_pos figure from R/48 (which used the wrong gating path); this pulls directly
# from cells/Singlets/CD45+/PBMC vivas/CD3+/PD-1+ in the same FlowJo export used for PD-1 MFI above.
pd1pc_raw <- pd1_raw %>% mutate(sample = str_extract(name, "^[^/]+"), path = str_remove(name, "^[^/]+/?"))
pd1pc <- pd1pc_raw %>% filter(path == "cells/Singlets/CD45+/PBMC vivas/CD3+/PD-1+") %>%
  filter(!str_detect(sample, "^SIN MARCAR|^Tumorales\\.fcs")) %>%
  transmute(sample, pct = as.numeric(stat), n = ncells) %>%
  rowwise() %>% mutate(parse_sample(sample)) %>% ungroup() %>%
  filter(condition != "Control") %>%
  mutate(condition = factor(condition, cond_lv) %>% droplevels(),
         environment = factor(environment, c("Basal", "Inflammatory")))

readr::write_csv(pd1pc, file.path("output/tidy", "pd1_pos_new.csv"))

pd1pc_summ <- pd1pc %>% group_by(environment, condition, time) %>%
  summarise(M_pct = mean(pct, na.rm = TRUE), SD_pct = sd(pct, na.rm = TRUE),
            M_n = mean(n, na.rm = TRUE), SD_n = sd(n, na.rm = TRUE), .groups = "drop")

panel2 <- function(df_pts, df_sum, y_pts, y_m, y_sd, ylab, title, offsets = cond_off2) {
  dd <- df_pts %>% mutate(x = xbase(environment, time) + offsets[as.character(condition)], val = .data[[y_pts]])
  ss <- df_sum %>% mutate(x = xbase(environment, time) + offsets[as.character(condition)],
                          M = .data[[y_m]], SD = .data[[y_sd]])
  ggplot() +
    geom_col(data = ss, aes(x, M, colour = condition), fill = NA, width = 0.28, linewidth = 1.0) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + replace_na(SD, 0), colour = condition),
                  width = 0.13, linewidth = 0.5) +
    geom_point(data = dd, aes(x, val, colour = condition), size = 1.6, alpha = 0.85) +
    geom_vline(xintercept = (xbase("Basal", 96) + xbase("Inflammatory", 24)) / 2,
              colour = "black", linewidth = 1.0) +
    scale_colour_manual(values = cond_pal, labels = cond_lab, name = NULL) +
    scale_x_continuous(breaks = c(xbase("Basal", c(24, 48, 96)), xbase("Inflammatory", c(24, 48, 96))),
                       labels = rep(c("24 h", "48 h", "96 h"), 2),
                       limits = c(tidx(24) - 0.6, xbase("Inflammatory", 96) + 0.6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = ylab, title = title) +
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

p_pct <- panel2(pd1pc, pd1pc_summ, "pct", "M_pct", "SD_pct", "% of CD3+ live PBMC", "PD-1+ (of CD3+ live PBMC)")
p_n   <- panel2(pd1pc, pd1pc_summ, "n",   "M_n",   "SD_n",   "Absolute count",      "PD-1+ (of CD3+ live PBMC)")

gl <- ggplotGrob(p_pct + theme(legend.position = "bottom", legend.direction = "horizontal"))
legend <- gl$grobs[[which(vapply(gl$grobs, function(x) x$name, "") == "guide-box")]]

png(file.path(OUT, "Fig_pd1_pos_new.png"), width = 11.5, height = 5.8, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit.c(unit(1, "null"), unit(2.2, "lines")))))
prt <- function(p, r, c) { pushViewport(viewport(layout.pos.row = r, layout.pos.col = c)); print(p, newpage = FALSE); popViewport() }
prt(p_pct + theme(legend.position = "none"), 1, 1)
prt(p_n   + theme(legend.position = "none"), 1, 2)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2)); grid.draw(legend); popViewport()
popViewport(); dev.off()
message("OK Fig_pd1_pos_new.png")
