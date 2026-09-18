# CD4+ and CD8+ T cells positive for HLA-DR (% of parent), from "CD4-CD8 DR+.xlsx".
# 4 independent panels in one file (rows = CD4 / CD8, cols = Basal / cytokine cocktail), same
# hollow-bar format as Fig6_activation_CD64CD11b: x = time (24/48/96 h), dodge = Resting (hollow) /
# Activated (solid fill), points + SD (n = 4 donors). Colour = population identity, matching
# Fig7C_unified (CD4 = #6a3d9a, CD8 = #cc4c02). Column code: CK = cytokine cocktail (Inflammatory),
# Beads = activation beads (Activated); no CK = Basal, no Beads = Resting.
source("R/fig_theme.R")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2); library(grid)
})

POP_COL <- c("CD4" = "#6a3d9a", "CD8" = "#cc4c02")
cond_lv  <- c("Resting", "Activated")
cond_off <- c(Resting = -0.19, Activated = 0.19)
env_lab  <- c(Basal = "Basal", Inflammatory = "TNF α, IL-α, IL-1 β")
tidx <- function(t) match(t, c(24, 48, 96))

raw <- read_excel("/mnt/c/Users/57319/Downloads/CD4-CD8 DR+.xlsx", sheet = "Hoja1", col_names = FALSE)
hdr <- as.character(raw[1, -1])
cmeta <- tibble(col = seq_along(hdr) + 1, label = hdr) %>%
  mutate(environment = ifelse(grepl("_CK_", label), "Inflammatory", "Basal"),
         condition = ifelse(grepl("Beads", label), "Activated", "Resting"),
         time = as.numeric(sub(".*_(\\d+)h$", "\\1", label)))

body <- raw[-1, ]
pop_of <- ifelse(grepl("^CD4", as.character(body[[1]])), "CD4", "CD8")
donor_of <- sub(".*_(\\d+)$", "D\\1", as.character(body[[1]]))

d0 <- bind_rows(lapply(seq_len(nrow(cmeta)), function(i) {
  cm <- cmeta[i, ]
  tibble(pop = pop_of, donor = donor_of, environment = cm$environment,
         condition = cm$condition, time = cm$time, value = as.numeric(body[[cm$col]]))
})) %>%
  mutate(condition = factor(condition, cond_lv),
         environment = factor(environment, c("Basal", "Inflammatory")),
         x = tidx(time) + cond_off[as.character(condition)])
s0 <- d0 %>% group_by(pop, environment, condition, time, x) %>%
  summarise(M = mean(value), SD = sd(value), .groups = "drop")

panel <- function(p, env) {
  dd <- filter(d0, pop == p, environment == env); ss <- filter(s0, pop == p, environment == env)
  col <- POP_COL[[p]]
  ggplot() +
    geom_col(data = ss, aes(x, M, alpha = condition), colour = col, fill = col, width = 0.34,
             linewidth = 1.1) +
    geom_errorbar(data = ss, aes(x, ymin = M, ymax = M + SD), colour = col, width = 0.14, linewidth = 0.55) +
    geom_point(data = dd, aes(x, value), colour = col, size = 1.3, alpha = 0.85,
               position = position_jitter(width = 0.05, height = 0)) +
    scale_alpha_manual(values = c(Resting = 0, Activated = 1), guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = c("24 h", "48 h", "96 h")) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Time (h)", y = paste0(p, "+ HLA-DR+ (%)"), title = env_lab[[as.character(env)]]) +
    theme_classic(base_size = 17) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.6),
          axis.ticks = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 14),
          axis.title = element_text(face = "bold", size = 17),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          panel.grid = element_blank())
}

p_cd4_b <- panel("CD4", "Basal"); p_cd4_i <- panel("CD4", "Inflammatory")
p_cd8_b <- panel("CD8", "Basal"); p_cd8_i <- panel("CD8", "Inflammatory")

# per-population legend: hollow = Resting, solid = Activated, in THAT population's own colour
make_legend <- function(p) {
  leg_df <- tibble(condition = factor(cond_lv, cond_lv), x = 1, y = 1)
  pleg <- ggplot(leg_df, aes(x, y, alpha = condition)) +
    geom_col(colour = POP_COL[[p]], fill = POP_COL[[p]], width = 0.6) +
    scale_alpha_manual(values = c(Resting = 0, Activated = 1), name = NULL,
                       labels = c("Resting PBMC", "Activated PBMC")) +
    theme_void(base_size = 15) + theme(legend.position = "right")
  gl <- ggplotGrob(pleg)
  gl$grobs[[which(vapply(gl$grobs, function(x) x$name, "") == "guide-box")]]
}
legend_cd4 <- make_legend("CD4"); legend_cd8 <- make_legend("CD8")

png("output/figuras/Fig_hladr_CD4CD8.png", width = 11.2, height = 10.0, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2,
             heights = unit.c(unit(1, "null"), unit(2, "lines"), unit(1, "null"), unit(2, "lines")))))
prt <- function(p, r, c) { pushViewport(viewport(layout.pos.row = r, layout.pos.col = c)); print(p, newpage = FALSE); popViewport() }
prt(p_cd4_b, 1, 1); prt(p_cd4_i, 1, 2)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2)); grid.draw(legend_cd4); popViewport()
prt(p_cd8_b, 3, 1); prt(p_cd8_i, 3, 2)
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1:2)); grid.draw(legend_cd8); popViewport()
popViewport(); dev.off()
message("OK Fig_hladr_CD4CD8.png")
