# R/15_fig_morfologia_barras.R — morphology bar charts, Prism style, split X axis.
# LEFT = Basal, RIGHT = Inflammatory; each with the 4 timepoints.
# Bars by PBMC condition; individual points; mean +/- SD; Tukey HSD brackets.
# No legend: each bar labelled (diagonal) with its condition; time as group header;
# environment titles on top; dashed divider; sci-notation Y for Area/Diameter.
source("R/fig_theme.R")

cond3_lv  <- c("Control (no PBMC)", "Resting PBMC", "Activated PBMC")
cond3_pal <- c("Control (no PBMC)" = "#9e9e9e", "Resting PBMC" = "#4292c6", "Activated PBMC" = "#e6550d")
cond_off  <- c("Control (no PBMC)" = -0.25, "Resting PBMC" = 0, "Activated PBMC" = 0.25)
cond_short <- c("Control (no PBMC)" = "Control", "Resting PBMC" = "PBMC resting", "Activated PBMC" = "PBMC activated")
stars <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", NA))))
sci_lab <- function(x) parse(text = ifelse(x == 0, "0",
                        gsub("e\\+?0*", "%*%10^", formatC(x, format = "e", digits = 1))))
tidx <- function(t) match(t, c(0, 24, 48, 96))
xbase <- function(env, t) ifelse(env == "Basal", tidx(t), tidx(t) + 5)  # Basal left, Inflammatory right

theme_prism <- function(base = 18) theme_classic(base_size = base) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = base - 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base - 4),
        axis.title = element_text(face = "bold", size = base + 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base + 5),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey35", size = base - 4),
        legend.position = "none", panel.grid = element_blank())

fig_morf_bar <- function(metric_key, ylab, file, title, sci = TRUE, bare = FALSE) {
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

  # Significance from the mixed-model post-hoc (emmeans Tukey; R/16). Stacked brackets.
  cmap <- c(None = "Control (no PBMC)", Resting = "Resting PBMC", Activated = "Activated PBMC")
  ph <- read_csv("output/estadistica/morfologia_posthoc.csv", show_col_types = FALSE) %>%
    filter(metric == metric_key, p.value < 0.05) %>%
    mutate(timef = as.numeric(as.character(timef)))
  br <- list()
  for (e in c("Basal","Inflammatory")) for (t in sort(unique(d$time))) {
    sub_ph <- ph %>% filter(environment == e, timef == t)
    if (nrow(sub_ph) == 0) next
    base <- max((s %>% filter(env == e, time == t) %>% mutate(top = M + SD))$top, na.rm = TRUE)
    sig <- list()
    for (i in seq_len(nrow(sub_ph))) {
      cs <- trimws(strsplit(sub_ph$contrast[i], "-")[[1]])
      x1 <- xbase(e,t) + cond_off[cmap[cs[1]]]; x2 <- xbase(e,t) + cond_off[cmap[cs[2]]]
      sig[[length(sig)+1]] <- data.frame(x = min(x1,x2), xend = max(x1,x2),
                                         lab = stars(sub_ph$p.value[i]), span = abs(x2 - x1))
    }
    sg <- bind_rows(sig) %>% arrange(span)
    sg$y <- base * (1 + 0.09 * seq_len(nrow(sg)))
    br[[length(br)+1]] <- sg
  }
  brdf <- if (length(br)) bind_rows(br) else NULL
  topbar <- max(s$M + s$SD, na.rm = TRUE)
  maxc <- if (!is.null(brdf)) max(brdf$y) else topbar
  time_y <- maxc * 1.11; env_y <- time_y * 1.07; ytop <- env_y * 1.05

  axis_df <- s %>% distinct(xpos, cond) %>% mutate(lab = cond_short[as.character(cond)]) %>% arrange(xpos)
  time_df <- s %>% distinct(env, time) %>% mutate(x = xbase(env, time), y = time_y, lab = paste0(time, " h"))

  p <- ggplot(s, aes(xpos, M, fill = cond)) +
    geom_col(width = 0.22, colour = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = M, ymax = M + SD), width = 0.1, linewidth = 0.4) +
    geom_point(data = d, aes(xpos, value), inherit.aes = FALSE,
               position = position_jitter(width = 0.05, height = 0), size = 0.9, colour = "black", alpha = 0.75) +
    geom_vline(xintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text(data = time_df, aes(x, y, label = lab), inherit.aes = FALSE, fontface = "bold", size = 5.4) +
    annotate("text", x = 2.5, y = env_y, label = "Basal", fontface = "bold", size = 7.6) +
    annotate("text", x = 7.5, y = env_y, label = "Inflammatory", fontface = "bold", size = 7.6) +
    scale_fill_manual(values = cond3_pal) +
    scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
    scale_y_continuous(labels = if (sci) sci_lab else waiver(),
                       expand = expansion(mult = c(0, 0.02)), limits = c(0, ytop)) +
    labs(title = if (bare) NULL else title,
         subtitle = if (bare) NULL else "Mean ± SD, n = 4; linear mixed model, Tukey-adjusted post-hoc (* p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001)",
         x = NULL, y = ylab) +
    theme_prism()
  if (!is.null(brdf)) {
    p <- p +
      geom_segment(data = brdf, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = brdf, aes(x = (x + xend)/2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.05, size = 5.6)
  }
  save_fig(p, file, w = 17, h = 6.6)
  message("OK ", file)
}

fig_morf_bar("Area", expression(bold(Area~(µm^2))), "Fig1_area_bars", "Spheroid area", sci = TRUE)
fig_morf_bar("Diametro", expression(bold(Diameter~(µm))), "Fig1_diameter_bars", "Spheroid diameter", sci = TRUE)
fig_morf_bar("Circularidad", "Circularity", "Fig1_circularity_bars", "Spheroid circularity", sci = FALSE)

# bare versions (no title/subtitle) for the unified vertical composition (A/B/C)
fig_morf_bar("Area", expression(bold(Area~(µm^2))), "Fig1c_area", "", sci = TRUE, bare = TRUE)
fig_morf_bar("Diametro", expression(bold(Diameter~(µm))), "Fig1c_diameter", "", sci = TRUE, bare = TRUE)
fig_morf_bar("Circularidad", "Circularity", "Fig1c_circularity", "", sci = FALSE, bare = TRUE)
