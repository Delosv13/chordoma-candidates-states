#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Test whether core biological-axis directions are robust to score scaling.
#
# Main inputs:
#   - data/processed/GSE230168_log2cpm_entrez.csv
#   - data/processed/GSE276715_log2tpm_entrez.csv
#   - data/processed/GSE216418_log2cpm_entrez.csv
#   - data/processed/GSE277707_log2norm_entrez.csv
#   - matching curated metadata files in data/processed/
#
# Main outputs:
#   - results/*_axis_z_scores_for_robustness.csv
#   - results/*_axis_rank_scores_for_robustness.csv
#   - results/signature_robustness_all_score_tests.csv
#   - results/signature_robustness_summary.csv
#   - results/signature_robustness_concordance_rate.csv
#
# Statistical approach:
#   Recomputes selected axes with both mean z-score and within-sample
#   percentile-rank scoring, then compares directions using Wilcoxon/BH-FDR.
#
# Interpretation notes:
#   This is a sensitivity analysis for directionality. Mixed concordance should
#   be reported as supplemental robustness evidence, not as complete replication.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

read_processed_matrix <- function(path) {
  x <- read_csv(path, show_col_types = FALSE)
  mat <- as.matrix(x[, -1, drop = FALSE])
  rownames(mat) <- as.character(x$entrez_id)
  storage.mode(mat) <- "numeric"
  mat
}

axis_gene_sets <- function() {
  sets <- selected_gene_sets()
  list(
    AXIS_IMMUNE_INTERFERON = unique(unlist(sets[c(
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_IL6_JAK_STAT3_SIGNALING",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "CUSTOM_ANTIGEN_PRESENTATION"
    )])),
    AXIS_JAK_STAT = sets[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]],
    AXIS_CELL_CYCLE = unique(unlist(sets[c(
      "HALLMARK_E2F_TARGETS",
      "HALLMARK_G2M_CHECKPOINT",
      "CUSTOM_CDK46_AXIS",
      "CUSTOM_PLK1_AXIS"
    )])),
    AXIS_CDK46 = sets[["CUSTOM_CDK46_AXIS"]],
    AXIS_PLK1 = sets[["CUSTOM_PLK1_AXIS"]],
    AXIS_RTK = sets[["CUSTOM_DRUGGABLE_RTKS"]],
    AXIS_CHORDOMA_IDENTITY = sets[["CUSTOM_CHORDOMA_MARKERS"]]
  )
}

axis_labels <- c(
  AXIS_IMMUNE_INTERFERON = "Immune/interferon",
  AXIS_JAK_STAT = "JAK/STAT",
  AXIS_CELL_CYCLE = "Cell-cycle",
  AXIS_CDK46 = "CDK4/6",
  AXIS_PLK1 = "PLK1",
  AXIS_RTK = "RTK",
  AXIS_CHORDOMA_IDENTITY = "Chordoma identity"
)

score_gene_sets_rank <- function(log_expr, gene_sets) {
  rank_expr <- apply(log_expr, 2, function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    (r - 0.5) / sum(!is.na(r))
  })
  rownames(rank_expr) <- rownames(log_expr)

  score_one <- function(genes) {
    present <- intersect(as.character(genes), rownames(rank_expr))
    if (length(present) < 3) {
      return(rep(NA_real_, ncol(rank_expr)))
    }
    colMeans(rank_expr[present, , drop = FALSE], na.rm = TRUE)
  }

  score_matrix <- do.call(rbind, lapply(gene_sets, score_one))
  colnames(score_matrix) <- colnames(rank_expr)
  out <- as.data.frame(t(score_matrix)) %>% rownames_to_column("sample")
  score_cols <- setdiff(colnames(out), "sample")
  out[score_cols] <- lapply(out[score_cols], function(x) as.numeric(scale(x)))
  out
}

run_dataset <- function(dataset, expr_path, meta_path, group_col, group_a, group_b) {
  expr <- read_processed_matrix(expr_path)
  meta <- read_csv(meta_path, show_col_types = FALSE)
  common <- intersect(colnames(expr), meta$sample)
  expr <- expr[, common, drop = FALSE]
  meta <- meta %>% filter(sample %in% common)
  sets <- axis_gene_sets()

  z_scores <- score_gene_sets(expr, sets)
  rank_scores <- score_gene_sets_rank(expr, sets)

  z_stats <- compare_scores(z_scores, meta, group_col, group_a, group_b) %>%
    transmute(
      dataset,
      contrast = paste0(group_a, "_vs_", group_b),
      axis = pathway,
      scoring = "z_score",
      n_a,
      n_b,
      mean_a,
      mean_b,
      delta = delta_a_minus_b,
      p_value,
      fdr
    )
  rank_stats <- compare_scores(rank_scores, meta, group_col, group_a, group_b) %>%
    transmute(
      dataset,
      contrast = paste0(group_a, "_vs_", group_b),
      axis = pathway,
      scoring = "rank_score",
      n_a,
      n_b,
      mean_a,
      mean_b,
      delta = delta_a_minus_b,
      p_value,
      fdr
    )

  write_csv(
    z_scores %>% left_join(meta, by = "sample") %>% mutate(dataset = dataset),
    paste0("results/", dataset, "_axis_z_scores_for_robustness.csv")
  )
  write_csv(
    rank_scores %>% left_join(meta, by = "sample") %>% mutate(dataset = dataset),
    paste0("results/", dataset, "_axis_rank_scores_for_robustness.csv")
  )

  bind_rows(z_stats, rank_stats)
}

message("Running signature robustness analyses")
stats <- bind_rows(
  run_dataset(
    "GSE230168",
    "data/processed/GSE230168_log2cpm_entrez.csv",
    "data/processed/GSE230168_metadata_curated.csv",
    "type",
    "chordoma_I",
    "chordoma_C"
  ),
  run_dataset(
    "GSE276715",
    "data/processed/GSE276715_log2tpm_entrez.csv",
    "data/processed/GSE276715_metadata_curated.csv",
    "condition",
    "T_DARPin",
    "control"
  ),
  run_dataset(
    "GSE216418",
    "data/processed/GSE216418_log2cpm_entrez.csv",
    "data/processed/GSE216418_metadata_curated.csv",
    "treatment",
    "palbociclib",
    "untreated"
  )
) %>%
  mutate(axis_label = axis_labels[axis]) %>%
  relocate(axis_label, .after = axis)

summary <- stats %>%
  select(dataset, contrast, axis, axis_label, scoring, delta, p_value, fdr, n_a, n_b) %>%
  pivot_wider(
    names_from = scoring,
    values_from = c(delta, p_value, fdr, n_a, n_b),
    names_glue = "{.value}_{scoring}"
  ) %>%
  mutate(
    direction_z_score = sign(delta_z_score),
    direction_rank_score = sign(delta_rank_score),
    direction_concordant = !is.na(direction_z_score) &
      !is.na(direction_rank_score) &
      direction_z_score == direction_rank_score,
    abs_delta_z_score = abs(delta_z_score),
    abs_delta_rank_score = abs(delta_rank_score)
  ) %>%
  arrange(dataset, axis)

primary_axes <- c(
  "AXIS_IMMUNE_INTERFERON",
  "AXIS_JAK_STAT",
  "AXIS_CELL_CYCLE",
  "AXIS_CDK46",
  "AXIS_PLK1",
  "AXIS_RTK",
  "AXIS_CHORDOMA_IDENTITY"
)
primary_summary <- summary %>%
  filter(axis %in% primary_axes) %>%
  summarise(
    n_tests = n(),
    n_concordant = sum(direction_concordant, na.rm = TRUE),
    concordance_rate = n_concordant / n_tests,
    include_in_main_text = concordance_rate >= 0.80,
    .groups = "drop"
  )

write_csv(stats, "results/signature_robustness_all_score_tests.csv")
write_csv(summary, "results/signature_robustness_summary.csv")
write_csv(primary_summary, "results/signature_robustness_concordance_rate.csv")

writeLines(c(
  "# Signature robustness analysis",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Method",
  "- Existing z-score gene-set scores were compared against rank-based gene-set scores.",
  "- The rank-based score ranks genes within each sample, averages percentile ranks across each axis, and z-scales the resulting axis score across samples.",
  "- Concordance is defined as matching effect direction for z-score and rank-score contrasts.",
  "",
  "## Main result",
  paste0(
    "- Primary-axis concordance: ",
    primary_summary$n_concordant,
    "/",
    primary_summary$n_tests,
    " (",
    round(100 * primary_summary$concordance_rate, 1),
    "%)."
  ),
  paste0("- Include robustness claim in main text: ", primary_summary$include_in_main_text, "."),
  "",
  "## Main outputs",
  "- results/signature_robustness_summary.csv",
  "- results/signature_robustness_concordance_rate.csv",
  "- results/*_axis_rank_scores_for_robustness.csv"
), "results/signature_robustness_summary.md")

message("Done: signature robustness analysis")
