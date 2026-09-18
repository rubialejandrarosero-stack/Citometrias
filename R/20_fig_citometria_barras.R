# R/20_fig_citometria_barras.R — flow-cytometry figures (Prism bar style, as morphology):
# absolute counts (main populations) + percentages (subpopulations & activation markers),
# split X by environment, time groups, PBMC resting/activated bars, individual points,
# mean +/- SD, significance from a linear mixed model (env x activation x time + (1|donor)).
# Also saves the full statistics (Type III ANOVA + Tukey post-hoc).
source("R/fig_theme.R")
suppressPackageStartupMessages({library(lmerTest); library(emmeans); library(scales)})
dir.create("output/estadistica", showWarnings = FALSE, recursive = TRUE)

cond2_lv  <- c("Resting", "Activated")
cond2_pal <- c("Resting" = "#4292c6", "Activated" = "#e6550d")
cond_off  <- c("Resting" = -0.19, "Activated" = 0.19)
cond_short <- c("Resting" = "PBMC resting", "Activated" = "PBMC activated")
stars <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", NA))))
tidx  <- function(t) match(t, c(24, 48, 96))
xbase <- function(env, t) ifelse(env == "Basal", tidx(t), tidx(t) + 4)

theme_prism <- function(base = 18) theme_classic(base_size = base) +
  theme(axis.line = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = base - 2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = base - 5),
        axis.title = element_text(face = "bold", size = base + 1),
        plot.title = element_text(face = "bold", hjust = 0.5, size = base + 3),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey35", size = base - 5),
        legend.position = "none", panel.grid = element_blank())

d0 <- read_csv("output/Infiltracion/infiltracion_conteos_porcentajes.csv", show_col_types = FALSE) %>%
  mutate(environment = factor(ifelse(Condicion == "INF", "Inflammatory", "Basal"),
                              levels = c("Basal", "Inflammatory")),
         cond = factor(ifelse(Activacion == "ACT", "Activated", "Resting"), levels = cond2_lv),
         timef = factor(Tiempo), donor = factor(Donante), time = Tiempo)

anova_txt <- character(0); posthoc_all <- list()

fig_cyto_bar <- function(pop, ylab, title, file, logit = FALSE, no_title = FALSE) {
  d <- d0 %>% mutate(value = .data[[pop]]) %>% filter(!is.na(value)) %>%
    mutate(mval = if (logit) qlogis(pmin(pmax(value / 100, 0.005), 0.995)) else value)
  res <- tryCatch({
    m <- lmer(mval ~ environment * cond * timef + (1 | donor), data = d)
    list(av = anova(m),
         ph = as.data.frame(pairs(emmeans(m, ~ cond | environment * timef), adjust = "tukey")))
  }, error = function(e) NULL)
  if (!is.null(res)) {
    anova_txt <<- c(anova_txt, paste0("\n================ ", title, "  (", pop, ") ================"),
                    utils::capture.output(print(res$av)))
    posthoc_all[[pop]] <<- res$ph %>% mutate(variable = pop, population = title)
  }

  s <- d %>% group_by(environment, time, cond) %>%
    summarise(M = mean(value), SD = sd(value), .groups = "drop") %>%
    mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])
  d <- d %>% mutate(xpos = xbase(environment, time) + cond_off[as.character(cond)])

  br <- NULL
  if (!is.null(res)) {
    ybase <- s %>% group_by(environment, time) %>%
      summarise(yb = max(M + SD, na.rm = TRUE), .groups = "drop") %>%
      mutate(environment = as.character(environment))
    br <- res$ph %>%
      mutate(environment = as.character(environment), time = as.numeric(as.character(timef)),
             lab = stars(p.value)) %>% filter(!is.na(lab)) %>%
      left_join(ybase, by = c("environment", "time")) %>%
      mutate(x = xbase(environment, time) + cond_off[["Resting"]],
             xend = xbase(environment, time) + cond_off[["Activated"]], y = yb * 1.06)
  }
  topbar <- max(s$M + s$SD, na.rm = TRUE)
  maxc <- if (!is.null(br) && nrow(br)) max(br$y) else topbar
  time_y <- maxc * 1.11; env_y <- time_y * 1.07; ytop <- env_y * 1.06

  axis_df <- s %>% distinct(xpos, cond) %>% mutate(lab = cond_short[as.character(cond)]) %>% arrange(xpos)
  time_df <- s %>% distinct(environment, time) %>% mutate(x = xbase(environment, time), y = time_y)

  p <- ggplot(s, aes(xpos, M, fill = cond)) +
    geom_col(width = 0.32, colour = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = M, ymax = M + SD), width = 0.14, linewidth = 0.4) +
    geom_point(data = d, aes(xpos, value), inherit.aes = FALSE,
               position = position_jitter(width = 0.06, height = 0), size = 1, colour = "black", alpha = 0.75) +
    geom_vline(xintercept = 4, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text(data = time_df, aes(x, y, label = time), inherit.aes = FALSE, fontface = "bold", size = 5) +
    annotate("text", x = 2, y = env_y, label = "Basal", fontface = "bold", size = 7) +
    annotate("text", x = 6, y = env_y, label = "Inflammatory", fontface = "bold", size = 7) +
    scale_fill_manual(values = cond2_pal) +
    scale_x_continuous(breaks = axis_df$xpos, labels = axis_df$lab) +
    scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.02)), limits = c(0, ytop)) +
    labs(title = if (no_title) NULL else title,
         subtitle = if (no_title) NULL else "Mean ± SD, n = 4; linear mixed model, Tukey post-hoc (* p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001)",
         x = NULL, y = ylab) +
    theme_prism()
  if (!is.null(br) && nrow(br)) {
    p <- p + geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = br, aes(x = (x + xend)/2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.05, size = 5)
  }
  save_fig(p, file, w = 12, h = 6)
  message("OK ", file)
}

# ---- absolute counts (main populations) ----
counts <- list(
  c("n_CD45",      "CD45+ total cell count", "Fig6c_CD45"),
  c("n_CD3",       "CD3+ cell count",        "Fig6c_CD3"),
  c("n_CD4",       "CD4+ cell count",        "Fig6c_CD4"),
  c("n_CD8",       "CD8+ cell count",        "Fig6c_CD8"),
  c("n_Mono_CD14", "CD14+ cell count",       "Fig6c_Monocytes"),
  c("n_B_CD19",    "CD19+ B cell count",     "Fig6c_Bcells"),
  c("n_NK_CD16",   "CD16+ NK cell count",    "Fig6c_NK"))
for (x in counts) fig_cyto_bar(x[1], x[2], NULL, x[3], no_title = TRUE)

# ---- percentages of CD45+ (subpopulations) ----
pcts <- list(
  c("pctCD45_CD3",  "CD3+ T cells",      "Fig6_CD3"),
  c("pctCD45_CD4",  "CD4+ T cells",      "Fig6_CD4"),
  c("pctCD45_CD8",  "CD8+ T cells",      "Fig6_CD8"),
  c("pctCD45_Mono", "Monocytes (CD14+)", "Fig6_Monocytes"))
for (x in pcts) fig_cyto_bar(x[1], "% of CD45+ leukocytes", x[2], x[3])
# rare subpopulations -> logit-transformed model (figure kept on raw % scale)
fig_cyto_bar("pctCD45_B",  "% of CD45+ leukocytes", "B cells (CD19+)",  "Fig6_Bcells", logit = TRUE)
fig_cyto_bar("pctCD45_NK", "% of CD45+ leukocytes", "NK cells (CD16+)", "Fig6_NK",     logit = TRUE)

# ---- activation markers (% of parent) ----
# HLA-DR-based activation removed (HLA-DR does not resolve in this panel; see R/33).
acts <- list(
  c("pctPar_CD64CD11b_de_Mono", "CD64+ CD11b+ (% of CD14+)",  "Fig6_act_CD64CD11b"))
for (x in acts) fig_cyto_bar(x[1], "% of parent population", x[2], x[3])

# ---- save statistics ----
writeLines(anova_txt, "output/estadistica/anova_citometria.txt")
bind_rows(posthoc_all) %>%
  select(population, variable, environment, timef, contrast, estimate, p.value) %>%
  mutate(p.value = round(p.value, 4)) %>%
  write_csv("output/estadistica/citometria_posthoc.csv")
message("Estadística guardada: output/estadistica/anova_citometria.txt + citometria_posthoc.csv")
