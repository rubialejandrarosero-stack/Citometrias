# R/11_fig_morfologia.R — Figure 1: spheroid morphology (area, diameter, circularity)
# Area & Diameter on log10 axis (scientific notation m x 10^n); Circularity linear.
source("R/fig_theme.R")
suppressPackageStartupMessages({library(gridExtra); library(grid); library(scales)})

d <- read_csv("output/tidy/morfologia_tidy.csv", show_col_types = FALSE) %>%
  mutate(group = make_group(environment, pbmc))

# drop non-physical area spikes (> 3x median)
d <- d %>% group_by(metric) %>%
  mutate(value = ifelse(metric == "Area" & value > 3 * median(value, na.rm = TRUE), NA, value)) %>%
  ungroup()

s <- d %>% filter(!is.na(value)) %>%
  group_by(metric, group, time) %>%
  summarise(M = mean(value), SE = sem(value), .groups = "drop")

# labels like  1.5 x 10^5
sci10 <- function(x) parse(text = gsub("e\\+?0*", " %*% 10^", formatC(x, format = "e", digits = 1)))

panel <- function(metric_name, ylab, log = FALSE, brks = waiver()) {
  p <- ggplot(filter(s, metric == metric_name), aes(time, M, colour = group, group = group)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 3, linewidth = 0.5) +
    scale_colour_manual(values = cond_pal, drop = FALSE) +
    scale_x_continuous(breaks = c(0, 24, 48, 96)) +
    labs(x = "Co-culture time (h)", y = ylab) +
    theme_frontiers() + theme(legend.position = "none")
  if (log) p <- p + scale_y_log10(breaks = brks, labels = sci10)
  p
}

pA <- panel("Area", "Area (a.u.)", log = TRUE, brks = c(9e4, 1.2e5, 1.5e5, 1.8e5))
pD <- panel("Diametro", "Diameter (µm)", log = TRUE, brks = c(350, 400, 450, 500))
pC <- panel("Circularidad", "Circularity")

# shared legend (extra spacing so long labels don't overlap)
lg_plot <- pC +
  theme(legend.position = "bottom",
        legend.text = element_text(margin = margin(r = 18, l = 3)),
        legend.key.width = unit(1.6, "lines"),
        legend.spacing.x = unit(0.5, "cm")) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
g <- ggplotGrob(lg_plot)
legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]

arr <- arrangeGrob(
  arrangeGrob(pA, pD, pC, nrow = 1),
  legend, nrow = 2, heights = c(10, 2),
  top = textGrob("Spheroid morphology over co-culture time",
                 gp = gpar(fontface = "bold", fontsize = 14), hjust = 0.5))

# wider panels (longer x axis) and shorter (shorter y axis)
dir.create("output/figuras", showWarnings = FALSE)
ggsave("output/figuras/Fig1_morphology.png", arr, width = 13, height = 4.0, dpi = 300)
ggsave("output/figuras/Fig1_morphology.pdf", arr, width = 13, height = 4.0)
message("OK Fig1_morphology (log axes on Area & Diameter)")
