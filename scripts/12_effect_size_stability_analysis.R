#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Quantify effect-size stability for prespecified core contrasts.
#
# Main inputs:
#   - results/GSE230168_pathway_scores.csv
#   - results/GSE230168_axis_scores.csv
#   - results/GSE276715_pathway_scores.csv
#   - results/GSE216418_pathway_scores.csv
#   - results/GSE216418_axis_scores.csv
#
# Main outputs:
#   - results/effect_size_stability_summary.csv
#   - results/bootstrap_pathway_delta_ci.csv
#   - results/leave_one_out_sensitivity.csv
#   - figures/Supplementary_Figure_effect_size_stability.*
#
# Statistical approach:
#   For each prespecified feature, computes mean deltas, Cliff's delta,
#   Wilcoxon p-values with BH-FDR, 1000 bootstrap confidence intervals, and
#   leave-one-out direction checks.
#
# Interpretation notes:
#   Stable core effects require direction consistency and bootstrap CIs that do
#   not cross zero; p-values alone are not treated as sufficient evidence.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

set.seed(20260608)

feature_labels <- c(
  HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6-JAK-STAT3",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon gamma",
  HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon alpha",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
  CUSTOM_DRUGGABLE_RTKS = "Druggable RTKs",
  CUSTOM_CHORDOMA_MARKERS = "Chordoma markers",
  HALLMARK_E2F_TARGETS = "E2F targets",
  HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
  CUSTOM_CDK46_AXIS = "CDK4/6 axis",
  CUSTOM_PLK1_AXIS = "PLK1 axis",
  immune_interferon = "Immune/interferon axis",
  cell_cycle = "Cell-cycle axis",
  chordoma_marker = "Chordoma marker axis"
)

analysis_plan <- tribble(
  ~dataset, ~score_path, ~group_col, ~case_value, ~control_value, ~feature,
  "GSE230168", "results/GSE230168_pathway_scores.csv", "type", "chordoma_I", "chordoma_C", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "GSE230168", "results/GSE230168_pathway_scores.csv", "type", "chordoma_I", "chordoma_C", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "GSE230168", "results/GSE230168_pathway_scores.csv", "type", "chordoma_I", "chordoma_C", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "GSE230168", "results/GSE230168_pathway_scores.csv", "type", "chordoma_I", "chordoma_C", "HALLMARK_INFLAMMATORY_RESPONSE",
  "GSE230168", "results/GSE230168_pathway_scores.csv", "type", "chordoma_I", "chordoma_C", "CUSTOM_DRUGGABLE_RTKS",
  "GSE230168", "results/GSE230168_axis_scores.csv", "type", "chordoma_I", "chordoma_C", "immune_interferon",
  "GSE230168", "results/GSE230168_axis_scores.csv", "type", "chordoma_I", "chordoma_C", "cell_cycle",
  "GSE230168", "results/GSE230168_axis_scores.csv", "type", "chordoma_I", "chordoma_C", "chordoma_marker",
  "GSE276715", "results/GSE276715_pathway_scores.csv", "condition", "T_DARPin", "control", "CUSTOM_CHORDOMA_MARKERS",
  "GSE276715", "results/GSE276715_pathway_scores.csv", "condition", "T_DARPin", "control", "HALLMARK_E2F_TARGETS",
  "GSE276715", "results/GSE276715_pathway_scores.csv", "condition", "T_DARPin", "control", "HALLMARK_G2M_CHECKPOINT",
  "GSE276715", "results/GSE276715_pathway_scores.csv", "condition", "T_DARPin", "control", "CUSTOM_CDK46_AXIS",
  "GSE276715", "results/GSE276715_pathway_scores.csv", "condition", "T_DARPin", "control", "CUSTOM_PLK1_AXIS",
  "GSE276715", "results/GSE276715_pathway_scores.csv", "condition", "T_DARPin", "control", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "GSE216418", "results/GSE216418_pathway_scores.csv", "treatment", "palbociclib", "untreated", "HALLMARK_E2F_TARGETS",
  "GSE216418", "results/GSE216418_pathway_scores.csv", "treatment", "palbociclib", "untreated", "HALLMARK_G2M_CHECKPOINT",
  "GSE216418", "results/GSE216418_pathway_scores.csv", "treatment", "palbociclib", "untreated", "CUSTOM_CDK46_AXIS",
  "GSE216418", "results/GSE216418_pathway_scores.csv", "treatment", "palbociclib", "untreated", "CUSTOM_PLK1_AXIS",
  "GSE216418", "results/GSE216418_axis_scores.csv", "treatment", "palbociclib", "untreated", "cell_cycle",
  "GSE216418", "results/GSE216418_axis_scores.csv", "treatment", "palbociclib", "untreated", "chordoma_marker"
)

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  comparisons <- outer(x, y, "-")
  mean(sign(comparisons))
}

summarise_vector <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  tibble(
    n_case = length(x),
    n_control = length(y),
    mean_case = mean(x),
    mean_control = mean(y),
    delta_case_minus_control = mean(x) - mean(y),
    cliffs_delta = cliffs_delta(x, y),
    wilcox_p = if (length(x) >= 2 && length(y) >= 2) {
      suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
    } else {
      NA_real_
    }
  )
}

bootstrap_ci <- function(x, y, n_boot = 1000) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(tibble(
      bootstrap_n = 0L,
      delta_ci_low = NA_real_,
      delta_ci_high = NA_real_,
      cliffs_ci_low = NA_real_,
      cliffs_ci_high = NA_real_
    ))
  }
  boot <- replicate(n_boot, {
    xb <- sample(x, length(x), replace = TRUE)
    yb <- sample(y, length(y), replace = TRUE)
    c(delta = mean(xb) - mean(yb), cliffs = cliffs_delta(xb, yb))
  })
  tibble(
    bootstrap_n = n_boot,
    delta_ci_low = unname(quantile(boot["delta", ], 0.025, na.rm = TRUE)),
    delta_ci_high = unname(quantile(boot["delta", ], 0.975, na.rm = TRUE)),
    cliffs_ci_low = unname(quantile(boot["cliffs", ], 0.025, na.rm = TRUE)),
    cliffs_ci_high = unname(quantile(boot["cliffs", ], 0.975, na.rm = TRUE))
  )
}

leave_one_out <- function(df, group_col, case_value, control_value, feature, full_delta) {
  bind_rows(lapply(df$sample, function(omitted_sample) {
    sub <- df %>% filter(sample != omitted_sample)
    x <- sub %>% filter(.data[[group_col]] == case_value) %>% pull(all_of(feature))
    y <- sub %>% filter(.data[[group_col]] == control_value) %>% pull(all_of(feature))
    loo <- summarise_vector(x, y)
    tibble(
      omitted_sample = omitted_sample,
      omitted_group = df[[group_col]][match(omitted_sample, df$sample)],
      loo_delta = loo$delta_case_minus_control,
      loo_cliffs_delta = loo$cliffs_delta,
      direction_flip = is.finite(full_delta) &&
        is.finite(loo$delta_case_minus_control) &&
        sign(full_delta) != 0 &&
        sign(loo$delta_case_minus_control) != sign(full_delta)
    )
  }))
}

run_one <- function(row) {
  score_path <- row$score_path
  feature <- row$feature
  if (!file.exists(score_path)) {
    stop("Missing score file: ", score_path)
  }
  df <- read_csv(score_path, show_col_types = FALSE)
  if (!feature %in% colnames(df)) {
    stop("Missing feature ", feature, " in ", score_path)
  }
  df <- df %>%
    filter(.data[[row$group_col]] %in% c(row$case_value, row$control_value)) %>%
    select(sample, all_of(row$group_col), all_of(feature)) %>%
    filter(!is.na(.data[[feature]]))

  x <- df %>% filter(.data[[row$group_col]] == row$case_value) %>% pull(all_of(feature))
  y <- df %>% filter(.data[[row$group_col]] == row$control_value) %>% pull(all_of(feature))
  base <- summarise_vector(x, y)
  ci <- bootstrap_ci(x, y)
  loo <- leave_one_out(df, row$group_col, row$case_value, row$control_value, feature, base$delta_case_minus_control) %>%
    mutate(
      dataset = row$dataset,
      contrast = paste0(row$case_value, "_vs_", row$control_value),
      feature = feature,
      feature_label = feature_labels[feature]
    ) %>%
    relocate(dataset, contrast, feature, feature_label)

  loo_summary <- loo %>%
    summarise(
      loo_n = n(),
      loo_delta_min = min(loo_delta, na.rm = TRUE),
      loo_delta_max = max(loo_delta, na.rm = TRUE),
      loo_direction_flips = sum(direction_flip, na.rm = TRUE),
      loo_direction_stable = loo_direction_flips == 0,
      .groups = "drop"
    )

  summary <- bind_cols(base, ci, loo_summary) %>%
    mutate(
      dataset = row$dataset,
      contrast = paste0(row$case_value, "_vs_", row$control_value),
      group_case = row$case_value,
      group_control = row$control_value,
      feature = feature,
      feature_label = feature_labels[feature],
      bootstrap_delta_excludes_zero = is.finite(delta_ci_low) &
        is.finite(delta_ci_high) &
        (delta_ci_low > 0 | delta_ci_high < 0),
      stable_core_effect = bootstrap_delta_excludes_zero & loo_direction_stable
    ) %>%
    relocate(dataset, contrast, group_case, group_control, feature, feature_label)

  list(summary = summary, loo = loo)
}

message("Running effect-size, bootstrap, and leave-one-out stability analyses")
results <- lapply(seq_len(nrow(analysis_plan)), function(i) run_one(analysis_plan[i, ]))
summary <- bind_rows(lapply(results, `[[`, "summary")) %>%
  mutate(fdr = p.adjust(wilcox_p, method = "BH")) %>%
  arrange(dataset, feature)
loo <- bind_rows(lapply(results, `[[`, "loo")) %>%
  arrange(dataset, feature, omitted_sample)

bootstrap_table <- summary %>%
  select(
    dataset, contrast, feature, feature_label, bootstrap_n,
    delta_case_minus_control, delta_ci_low, delta_ci_high,
    cliffs_delta, cliffs_ci_low, cliffs_ci_high,
    bootstrap_delta_excludes_zero
  )

write_csv(summary, "results/effect_size_stability_summary.csv")
write_csv(bootstrap_table, "results/bootstrap_pathway_delta_ci.csv")
write_csv(loo, "results/leave_one_out_sensitivity.csv")

plot_df <- summary %>%
  mutate(
    label = paste(dataset, feature_label, sep = ": "),
    label = factor(label, levels = rev(label)),
    stability = if_else(stable_core_effect, "Stable core effect", "Context-dependent/unstable")
  )

p_delta <- ggplot(plot_df, aes(x = delta_case_minus_control, y = label, color = stability)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_errorbar(aes(xmin = delta_ci_low, xmax = delta_ci_high), orientation = "y", width = 0.18, linewidth = 0.45) +
  geom_point(size = 2.1) +
  scale_color_manual(values = c("Stable core effect" = "#0072B2", "Context-dependent/unstable" = "#D55E00")) +
  labs(title = "a", x = "Mean score delta with 95% bootstrap CI", y = NULL, color = NULL) +
  theme_nature(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

p_cliff <- ggplot(plot_df, aes(x = cliffs_delta, y = label, color = stability)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_point(size = 2.1) +
  scale_color_manual(values = c("Stable core effect" = "#0072B2", "Context-dependent/unstable" = "#D55E00")) +
  labs(title = "b", x = "Cliff's delta", y = NULL, color = NULL) +
  theme_nature(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

fig <- plot_grid(p_delta, p_cliff, nrow = 1, rel_widths = c(1.25, 1))
ggsave("figures/Supplementary_Figure_effect_size_stability.png", fig, width = 13, height = 8.5, dpi = 300)
ggsave("figures/Supplementary_Figure_effect_size_stability.pdf", fig, width = 13, height = 8.5)

stable_count <- sum(summary$stable_core_effect, na.rm = TRUE)
writeLines(c(
  "# Effect-size and stability analysis",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Method",
  "- Main pathway and axis contrasts were summarized with mean score deltas, Cliff's delta, and Wilcoxon tests.",
  "- Bootstrap 95% confidence intervals used 1000 within-group resamples.",
  "- Leave-one-out sensitivity was used to flag comparisons whose effect direction changed after omitting one sample.",
  "",
  "## Main result",
  paste0("- Stable core effects: ", stable_count, "/", nrow(summary), "."),
  "- Stable core effects require a bootstrap delta CI excluding zero and no leave-one-out direction flips.",
  "",
  "## Main outputs",
  "- results/effect_size_stability_summary.csv",
  "- results/bootstrap_pathway_delta_ci.csv",
  "- results/leave_one_out_sensitivity.csv",
  "- figures/Supplementary_Figure_effect_size_stability.png"
), "results/effect_size_stability_summary.md")

message("Done: effect-size stability analysis")
