# R/31_fig_morfologia_dotplot.R — morphology as scatter dot plots (individual spheroids
# + mean ± SD), 3 conditions per timepoint, split X (Basal | Inflammatory). Y-axis zoomed
# to the data (not forced to 0) so differences are visible. No title; significance from
# the mixed-model post-hoc (R/16). For the unified vertical A/B/C composition.
source("R/fig_theme.R")
suppressPackageStartupMessages(library(scales))

cond3_lv  <- c("Control (no PBMC)", "Resting PBMC", "Activated PBMC")
cond3_pal <- c("Control (no PBMC)" = "#9e9e9e", "Resting PBMC" = "#4292c6", "Activated PBMC" = "#e6550d")
cond_off  <- c("Control (no PBMC)" = -0.25, "Resting PBMC" = 0, "Activated PBMC" = 0.25)
cond_short <- c("Control (no PBMC)" = "Control", "Resting PBMC" = "PBMC resting", "Activated PBMC" = "PBMC activated")
stars <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", NA))))
sci_lab <- function(x) parse(text = ifelse(x == 0, "0",
                        gsub("e\\+?0*", "%*%10^", formatC(x, format = "e", digits = 1))))
tidx <- function(t) match(t, c(0, 24, 48, 96))
xbase <- function(env, t) ifelse(env == "Basal", tidx(t), tidx(t) + 5)

theme_prism <- function(base = 21) theme_classic(base_size = base) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = base - 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base - 4),
        axis.title.y = element_text(face = "bold", size = base + 2),
        legend.position = "none", panel.grid = element_blank())

fig_morf_dot <- function(metric_key, ylab, file, sci = TRUE) {
  d <- read_csv("output/tidy/morfologia_tidy.csv", show_col_types = FALSE) %>%
    filter(metric == metric_key) %>%
    mutate(env = ifelse(environment == "Inflammatory", "Inflammatory", "Basal"),
           cond = factor(recode(pbmc, None = "Control (no PBMC)", Resting = "Resting PBMC",
                                Activated = "Activated PBMC"), levels = cond3_lv))
  if (metric_key == "Area")
    d <- d %>% mutate(value = ifelse(value > 3 * median(value, na.rm = TRUE), NA, value))
  d <- d %>% filter(!is.na(value)) %>% mutate(xpos = xbase(env, time) + cond_off[as.character(cond)])
  s <- d %>% group_by(env, time, cond) %>%
    summarise(M = mean(value), SD = sd(value), .groups = "drop") %>%
    mutate(xpos = xbase(env, time) + cond_off[as.character(cond)])

  cmap <- c(None = "Control (no PBMC)", Resting = "Resting PBMC", Activated = "Activated PBMC")
  ph <- read_csv("output/estadistica/morfologia_posthoc.csv", show_col_types = FALSE) %>%
    filter(metric == metric_key, p.value < 0.05) %>% mutate(timef = as.numeric(as.character(timef)))
  br <- list()
  for (e in c("Basal","Inflammatory")) for (t in sort(unique(d$time))) {
    sub_ph <- ph %>% filter(environment == e, timef == t)
    if (nrow(sub_ph) == 0) next
    base <- max((d %>% filter(env == e, time == t))$value, na.rm = TRUE)
    sig <- list()
    for (i in seq_len(nrow(sub_ph))) {
      cs <- trimws(strsplit(sub_ph$contrast[i], "-")[[1]])
      x1 <- xbase(e,t) + cond_off[cmap[cs[1]]]; x2 <- xbase(e,t) + cond_off[cmap[cs[2]]]
      sig[[length(sig)+1]] <- data.frame(x = min(x1,x2), xend = max(x1,x2),
                                         lab = stars(sub_ph$p.value[i]), span = abs(x2 - x1))
    }
    sg <- bind_rows(sig) %>% arrange(span)
    sg$y <- base * (1 + 0.085 * seq_len(nrow(sg)))
    br[[length(br)+1]] <- sg
  }
  brdf <- if (length(br)) bind_rows(br) else NULL

  ymin <- min(d$value, na.rm = TRUE); ymax_pt <- max(d$value, na.rm = TRUE)
  maxc <- if (!is.null(brdf)) max(brdf$y) else ymax_pt
  time_y <- maxc * 1.14; env_y <- time_y * 1.12; ytop <- env_y * 1.05
  ylo <- if (sci) ymin * 0.88 else 0  # circularity keeps 0 baseline

  axis_df <- s %>% distinct(xpos, cond) %>% mutate(lab = cond_short[as.character(cond)]) %>% arrange(xpos)
  time_df <- s %>% distinct(env, time) %>% mutate(x = xbase(env, time), y = time_y, lab = paste0(time, " h"))

  p <- ggplot() +
    geom_point(data = d, aes(xpos, value, colour = cond), shape = 1, size = 1.9, stroke = 0.9,
               position = position_jitter(width = 0.05, height = 0)) +
    geom_errorbar(data = s, aes(xpos, ymin = M - SD, ymax = M + SD, colour = cond), width = 0.12, linewidth = 0.6) +
    geom_errorbar(data = s, aes(xpos, ymin = M, ymax = M, colour = cond), width = 0.22, linewidth = 0.9) +
    geom_vline(xintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text(data = time_df, aes(x, y, label = lab), inherit.aes = FALSE, fontface = "bold", size = 6.2) +
    annotate("text", x = 2.5, y = env_y, label = "Basal", fontface = "bold", size = 8.4) +
    annotate("text", x = 7.5, y = env_y, label = "Inflammatory", fontface = "bold", size = 8.4) +
    scale_colour_manual(values = cond3_pal) +
    scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
    scale_y_continuous(labels = if (sci) sci_lab else waiver(),
                       expand = expansion(mult = c(0.02, 0.02))) +
    coord_cartesian(ylim = c(ylo, ytop)) +
    labs(x = NULL, y = ylab) +
    theme_prism()
  if (!is.null(brdf)) {
    p <- p + geom_segment(data = brdf, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = brdf, aes(x = (x + xend)/2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.05, size = 6.2)
  }
  save_fig(p, file, w = 17, h = 6.6)
  message("OK ", file)
}

fig_morf_dot("Area", expression(bold(Area~(µm^2))), "Fig1d_area", sci = TRUE)
fig_morf_dot("Diametro", expression(bold(Diameter~(µm))), "Fig1d_diameter", sci = TRUE)
fig_morf_dot("Circularidad", "Circularity", "Fig1d_circularity", sci = FALSE)
