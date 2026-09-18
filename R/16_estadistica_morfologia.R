# R/16_estadistica_morfologia.R
# Linear mixed-effects model for spheroid morphology.
#   value ~ environment * PBMC * time  + (1 | donor)
# donor = experimental replicate/batch (random effect; controls get batch 1-4).
# Type III ANOVA (Satterthwaite df) + Tukey-adjusted post-hoc (emmeans).
suppressPackageStartupMessages({
  library(lmerTest); library(emmeans); library(dplyr); library(readr); library(tidyr)
})
dir.create("output/estadistica", showWarnings = FALSE, recursive = TRUE)

d <- read_csv("output/tidy/morfologia_tidy.csv", show_col_types = FALSE) %>%
  group_by(metric, environment, pbmc, time) %>%
  mutate(donor = ifelse(is.na(donor), row_number(), donor)) %>% ungroup() %>%
  mutate(donor = factor(donor),
         environment = factor(environment, levels = c("Basal", "Inflammatory")),
         pbmc = factor(pbmc, levels = c("None", "Resting", "Activated")),
         timef = factor(time))

# drop non-physical area spike (as in the figure)
d <- d %>% group_by(metric) %>%
  mutate(value = ifelse(metric == "Area" & value > 3 * median(value, na.rm = TRUE), NA, value)) %>%
  ungroup() %>% filter(!is.na(value))

emmeans::emm_options(lmerTest.limit = 200, pbkrtest.limit = 200)
posthoc_all <- list(); sink("output/estadistica/anova_morfologia.txt")
for (mk in c("Area", "Diametro", "Circularidad")) {
  dd <- filter(d, metric == mk)
  m <- lmer(value ~ environment * pbmc * timef + (1 | donor), data = dd)
  cat("\n=====================================================================\n")
  cat("METRIC:", mk, "   (linear mixed model, random intercept = donor)\n")
  cat("=====================================================================\n")
  cat("\n-- Type III ANOVA (Satterthwaite) --\n"); print(anova(m))
  # post-hoc: PBMC conditions within each environment x time
  emm <- emmeans(m, ~ pbmc | environment * timef)
  pw <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    mutate(metric = mk) %>% select(metric, environment, timef, contrast, estimate, p.value)
  posthoc_all[[mk]] <- pw
}
sink()

posthoc <- bind_rows(posthoc_all)
write_csv(posthoc, "output/estadistica/morfologia_posthoc.csv")

# console summary focused on the question (Area) + significant post-hoc counts
cat("==== ANOVA saved to output/estadistica/anova_morfologia.txt ====\n\n")
cat("Significant post-hoc pairs (p<0.05) per metric:\n")
print(posthoc %>% filter(p.value < 0.05) %>% count(metric, name = "n_sig_pairs"))
cat("\nAREA — all PBMC pairwise p-values (should be non-significant):\n")
print(posthoc %>% filter(metric == "Area") %>%
        mutate(p.value = round(p.value, 3)) %>%
        arrange(p.value) %>% head(8))
