# CD14+CD64+CD11b+HLA-DR+ monocytes (% of parent), from "CD14+CD64+CD11b+HLA-DR.xlsx".
# 2 independent panels (Basal / cytokine cocktail), same hollow-bar format as Fig7F_hladr_CD4CD8:
# x = time (24/48/96 h), dodge = Resting (hollow) / Activated (solid fill), points + SD (n = 4
# donors). Colour: single population -> dark teal (distinct from CD4 purple / CD8 orange already
# used). Column code: CK = cytokine cocktail (Inflammatory), Beads = activation beads (Activated);
# no CK = Basal, no Beads = Resting.
source("R/fig_theme.R")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(grid)
})

cond_lv  <- c("Resting", "Activated")
cond_pal <- c(Resting = "#4292c6", Activated = "#e31a1c")   # blue / red, same as the rest of the family
cond_lab <- c(Resting = "Resting PBMC", Activated = "Activated PBMC")
cond_off <- c(Resting = -0.19, Activated = 0.19)
env_lab  <- c(Basal = "Basal", Inflammatory = "TNF α, IL-α, IL-1 β")
tidx <- function(t) match(t, c(24, 48, 96))

raw <- read_excel("/mnt/c/Users/57319/Downloads/CD14+CD64+CD11b+HLA-DR.xlsx", sheet = "Hoja1", col_names = FALSE)
hdr <- as.character(raw[1, -1])
cmeta <- tibble(col = seq_along(hdr) + 1, label = hdr) %>%
  mutate(environment = ifelse(grepl("_CK_", label), "Inflammatory", "Basal"),
         condition = ifelse(grepl("Beads", label), "Activated", "Resting"),
         time = as.numeric(sub(".*_(\\d+)h$", "\\1", label)))

body <- raw[-1, ]
donor_of <- sub(".*_(\\d+)$", "D\\1", as.character(body[[1]]))

d0 <- bind_rows(lapply(seq_len(nrow(cmeta)), function(i) {
  cm <- cmeta[i, ]
  tibble(donor = donor_of, environment = cm$environment, condition = cm$condition, time = cm$time,
         value = as.numeric(body[[cm$col]]))
})) %>%
  mutate(condition = factor(condition, cond_lv),
         environment = factor(environment, c("Basal", "Inflammatory")),
         x = tidx(time) + cond_off[as.character(condition)])
s0 <- d0 %>% group_by(environment, condition, time, x) %>%
  summarise(M = mean(value), SD = sd(value), .groups = "drop")

panel <- function(env) {
  dd <- filter(d0, environment == env); ss <- filter(s0, environment == env)
  ggplot() +
    geom_col(data = ss, aes(x, M, colour = condition), fill = NA, width = 0.34, linewidth = 1.1) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + SD, colour = condition), width = 0.14, linewidth = 0.55) +
    geom_point(data = dd, aes(x, value, colour = condition), size = 1.3, alpha = 0.85,
               position = position_jitter(width = 0.05, height = 0)) +
    scale_colour_manual(values = cond_pal, labels = cond_lab, name = NULL) +
    scale_x_continuous(breaks = 1:3, labels = c("24", "48", "96")) +
    scale_y_continuous(limits = c(0, 60), expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Time (h)", y = "CD14+CD64+CD11b+HLA-DR+ (%)", title = env_lab[[as.character(env)]]) +
    theme_classic(base_size = 17) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 14),
          axis.title = element_text(face = "bold", size = 17),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          legend.position = "none", panel.grid = element_blank())
}

p_b <- panel("Basal"); p_i <- panel("Inflammatory")

gl <- ggplotGrob(p_b + theme(legend.position = "bottom", legend.direction = "horizontal",
                             legend.text = element_text(size = 15)))
legend <- gl$grobs[[which(vapply(gl$grobs, function(x) x$name, "") == "guide-box")]]

png("output/figuras/Fig_hladr_monocytes.png", width = 11.2, height = 5.4, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit.c(unit(1, "null"), unit(2, "lines")))))
prt <- function(p, r, c) { pushViewport(viewport(layout.pos.row = r, layout.pos.col = c)); print(p, newpage = FALSE); popViewport() }
prt(p_b, 1, 1); prt(p_i, 1, 2)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2)); grid.draw(legend); popViewport()
popViewport(); dev.off()
message("OK Fig_hladr_monocytes.png")
