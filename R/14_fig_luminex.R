# R/14_fig_luminex.R — Figure 8: secretome (Luminex) heatmap + PCA
source("R/fig_theme.R")
suppressPackageStartupMessages({library(gridExtra); library(grid); library(tibble)})

d <- read_csv("output/tidy/luminex_tidy.csv", show_col_types = FALSE) %>%
  mutate(logv = log10(value + 1),
         pbmc_s = recode(pbmc, None = "Ctrl", Resting = "Rest", Activated = "Act"),
         col = factor(paste0(pbmc_s, " ", time, "h"),
                      levels = c(outer(c("Ctrl","Rest","Act"), c(24,48,96),
                                       function(a,b) paste0(a," ",b,"h")))))

# ---------- HEATMAP (z-score per cytokine, clustered rows, facet by environment) ----------
agg <- d %>% group_by(cytokine, environment, col) %>%
  summarise(m = mean(logv, na.rm = TRUE), .groups = "drop") %>%
  group_by(cytokine) %>% mutate(z = as.numeric(scale(m))) %>% ungroup()

mat <- agg %>% mutate(k = paste(environment, col)) %>%
  select(cytokine, k, z) %>% pivot_wider(names_from = k, values_from = z) %>%
  column_to_rownames("cytokine")
ord <- rownames(mat)[hclust(dist(as.matrix(mat)))$order]
agg$cytokine <- factor(agg$cytokine, levels = ord)

hm <- ggplot(agg, aes(col, cytokine, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  facet_wrap(~environment, nrow = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       name = "z-score") +
  labs(title = "Secreted cytokine/chemokine profile of tumor spheroids",
       subtitle = "Luminex multiplex, log10 mean per condition (row z-score, n = 4 donors)",
       x = NULL, y = NULL) +
  theme_frontiers(base = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right",
        legend.title = element_text(face = "bold"), panel.grid = element_blank())
save_fig(hm, "Fig8A_secretome_heatmap", w = 10, h = 8)
message("OK Fig8A_secretome_heatmap")

# ---------- PCA per environment ----------
wide <- d %>% mutate(donor = ifelse(is.na(donor), 0, donor),
                     sid = paste(environment, pbmc, time, donor, sep = "|")) %>%
  select(sid, environment, pbmc, time, cytokine, logv) %>%
  pivot_wider(names_from = cytokine, values_from = logv, values_fn = mean)

pca_panel <- function(envlab) {
  w <- wide %>% filter(environment == envlab)
  X <- w %>% select(-sid, -environment, -pbmc, -time) %>% as.matrix()
  X <- X[, apply(X, 2, function(c) sd(c, na.rm = TRUE) > 0)]
  X[is.na(X)] <- 0
  pc <- prcomp(X, scale. = TRUE)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2],
                   group = make_group(w$environment, w$pbmc), time = factor(w$time))
  ggplot(df, aes(PC1, PC2, colour = group, shape = time)) +
    geom_point(size = 2.6, alpha = 0.9) +
    scale_colour_manual(values = cond_pal, drop = TRUE) +
    labs(title = envlab, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
         shape = "Time (h)") +
    guides(colour = guide_legend(nrow = 3, order = 1), shape = guide_legend(nrow = 3, order = 2)) +
    theme_frontiers(base = 10) +
    theme(legend.position = "bottom", legend.title = element_text(face = "bold"),
          legend.box = "horizontal")
}
pB <- pca_panel("Basal"); pI <- pca_panel("Inflammatory")
arr <- arrangeGrob(pB, pI, nrow = 1,
                   top = textGrob("Principal component analysis of the spheroid secretome",
                                  gp = gpar(fontface = "bold", fontsize = 13)))
ggsave("output/figuras/Fig8B_secretome_PCA.png", arr, width = 12, height = 5.6, dpi = 300)
ggsave("output/figuras/Fig8B_secretome_PCA.pdf", arr, width = 12, height = 5.6)
message("OK Fig8B_secretome_PCA")
