# R/44_fig_cytokines_rest.R — remaining 22 cytokines from the Luminex secretome (pg/mL), i.e. all
# analytes NOT already shown in Fig_cytokines_key.png (R/43). Same design: independent panels per
# cytokine, 3 timepoints (24/48/96 h), Basal | Inflammatory groups (thick divider + in-panel
# "No inflamatorio"/"Inflamatorio" label), Control/Resting/Activated as dodged hollow bars + donor
# points + SD.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readr); library(scales); library(grid); library(tidyr)
})
OUT <- "output/figuras"
cond_lv  <- c("Control", "Resting", "Activated")
cond_pal <- c(Control = "black", Resting = "#4292c6", Activated = "#e31a1c")
cond_lab <- c("Control", "Resting PBMC", "Activated PBMC")
cond_off <- c(Control = -0.30, Resting = 0, Activated = 0.30)
tidx <- function(t) match(t, c(24, 48, 96))
xbase <- function(env, t) ifelse(env == "Basal", tidx(t), tidx(t) + 4)

KEY_DONE <- c("CD25", "CXCL10/IP10/CRG2", "CXCL-9", "IFN-gamma", "IL-10", "IL-17",
             "FGF basic/ FGF2/bFGF", "GM-CSF")
LAB <- c("CCL5 RANTES" = "CCL5", "CXCL-9" = "CXCL9", "CXCL10/IP10/CRG2" = "CXCL10",
        "IL-12/IL-23p40" = "IL-12/23", "FGF basic/ FGF2/bFGF" = "FGF-2",
        "IL-1beta" = "IL-1β", "TNF-alpha" = "TNF-α", "IFN-gamma" = "IFN-γ", "IFN-alpha" = "IFN-α")

all_raw <- read_csv("output/tidy/luminex_tidy.csv", show_col_types = FALSE) %>% distinct(cytokine) %>% pull(cytokine)
REST_RAW <- setdiff(all_raw, KEY_DONE)
KEY <- setNames(ifelse(REST_RAW %in% names(LAB), LAB[REST_RAW], REST_RAW), REST_RAW)
message("Plotting ", length(KEY), " remaining cytokines: ", paste(unname(KEY), collapse = ", "))

d0 <- read_csv("output/tidy/luminex_tidy.csv", show_col_types = FALSE) %>%
  filter(cytokine %in% names(KEY)) %>%
  mutate(condition = factor(ifelse(pbmc %in% c("None", NA), "Control", pbmc), cond_lv),
         cyt = factor(KEY[cytokine], levels = unname(KEY)),
         environment = factor(environment, c("Basal", "Inflammatory")),
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
    geom_col(data = ss, aes(x, M, colour = condition), fill = NA, width = 0.26, linewidth = 1.8) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + SD, colour = condition),
                  width = 0.15, linewidth = 1.1) +
    geom_point(data = dd, aes(x, value, colour = condition), size = 1.1, alpha = 0.8,
               position = position_jitter(width = 0.05, height = 0)) +
    geom_vline(xintercept = 4, colour = "black", linewidth = 1.0) +
    geom_text(data = env_labs, aes(x, y, label = lab), inherit.aes = FALSE,
              size = 12, fontface = "italic", colour = "grey30", vjust = 1) +
    scale_colour_manual(values = cond_pal, labels = cond_lab, name = NULL) +
    scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7), labels = rep(c("24", "48", "96"), 2),
                       limits = c(0.0, 7.85)) +
    scale_y_continuous(labels = label_comma(), limits = c(0, ytop),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Time (h)", y = "Concentration (pg/mL)", title = as.character(cy)) +
    theme_classic(base_size = 40) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.x = element_blank(),
          axis.text = element_text(colour = "black", size = 40),
          axis.text.x = element_text(size = 34, angle = 0, hjust = 0.5),
          axis.title.y = element_text(face = "bold", size = 40),
          axis.title.x = element_text(face = "bold", size = 40),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 40),
          legend.position = "top", legend.direction = "horizontal",
          legend.text = element_text(size = 40), legend.key.size = unit(1.8, "lines"),
          legend.margin = margin(1, 1, 1, 1), panel.grid = element_blank())

  p
}

cys <- unname(KEY)

# Cairo has a hard cap on raster surface size, and one 22-panel canvas at this font size (40pt,
# tuned for a 2-col x 4-row grid) would exceed it -- split into files of at most 4 panel-rows
# (8 cytokines) each, same 2-col layout/size as Fig_cytokines_key.png.
MAX_PER_FILE <- 8
n_files <- ceiling(length(cys) / MAX_PER_FILE)

render_file <- function(cys_chunk, outfile) {
  n_rows_panels <- ceiling(length(cys_chunk) / 2)
  heights_u <- rep(unit.c(unit(1, "null"), unit(2.4, "lines")), n_rows_panels)
  heights_u <- heights_u[seq_len(length(heights_u) - 1)]

  png(file.path(OUT, outfile), width = 34.0, height = 42.0 * n_rows_panels / 4,
     units = "in", res = 300, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(heights_u), 2, heights = heights_u)))
  prt <- function(p, r, c) {
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = c, clip = "off"))
    print(p, newpage = FALSE); popViewport()
  }
  panel_rows <- seq(1, by = 2, length.out = n_rows_panels)
  for (i in seq_along(cys_chunk)) {
    r <- panel_rows[(i - 1) %/% 2 + 1]
    prt(panel(cys_chunk[i]), r, (i - 1) %% 2 + 1)
  }
  popViewport(); dev.off()
  message("OK ", outfile)
}

for (f in seq_len(n_files)) {
  idx <- ((f - 1) * MAX_PER_FILE + 1):min(f * MAX_PER_FILE, length(cys))
  outfile <- if (n_files == 1) "Fig_cytokines_rest.png" else sprintf("Fig_cytokines_rest_%d.png", f)
  render_file(cys[idx], outfile)
}

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
