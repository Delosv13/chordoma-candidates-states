#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Add GSE230168 microenvironment-marker and biomarker context.
#
# Main inputs:
#   - data/processed/GSE230168_log2cpm_entrez.csv
#   - data/processed/GSE230168_metadata_curated.csv
#   - results/GSE230168_axis_scores.csv
#
# Main outputs:
#   - results/GSE230168_microenvironment_scores.csv
#   - results/GSE230168_microenvironment_score_tests.csv
#   - results/GSE230168_biomarker_panel.csv
#   - results/GSE230168_biomarker_group_tests.csv
#   - results/GSE230168_biomarker_axis_correlations.csv
#   - figures/Supplementary_Figure_GSE230168_biomarkers_microenvironment.*
#
# Statistical approach:
#   Bulk marker-set scores are computed by mean z-score. Group contrasts use
#   Wilcoxon/BH-FDR, and biomarker-axis associations use Spearman correlation.
#
# Interpretation notes:
#   Microenvironment scores are bulk RNA marker summaries, not cell-type
#   deconvolution.
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

microenvironment_sets <- function() {
  base_sets <- selected_gene_sets()
  list(
    MICRO_ANTIGEN_PRESENTATION = base_sets[["CUSTOM_ANTIGEN_PRESENTATION"]],
    MICRO_T_CELL_INFLAMMATION = c("915", "916", "914", "925", "926", "3002", "5551", "4818", "4283", "3627", "3458", "6772", "5133", "3902"),
    MICRO_MACROPHAGE_MYELOID = c("968", "9332", "1436", "7940", "2214", "3684", "712", "713", "714", "4481", "7305"),
    MICRO_CAF_STROMAL = c("59", "1277", "1278", "1281", "1289", "1634", "4060", "2191", "5159", "6876"),
    MICRO_ENDOTHELIAL_ANGIOGENIC = c("5175", "7450", "3791", "2321", "2022", "1003", "7422", "285"),
    MICRO_CHECKPOINT_CONTEXT = c("29126", "80380", "5133", "1493", "3902", "84868", "201633", "961", "3620", "64115")
  )
}

micro_labels <- c(
  MICRO_ANTIGEN_PRESENTATION = "Antigen presentation",
  MICRO_T_CELL_INFLAMMATION = "T-cell inflammation",
  MICRO_MACROPHAGE_MYELOID = "Macrophage/myeloid",
  MICRO_CAF_STROMAL = "CAF/stromal",
  MICRO_ENDOTHELIAL_ANGIOGENIC = "Endothelial/angiogenic",
  MICRO_CHECKPOINT_CONTEXT = "Checkpoint context"
)

biomarker_panel <- tribble(
  ~axis, ~gene_symbol, ~entrez_id,
  "Chordoma identity", "TBXT", "6862",
  "Cell-cycle", "MKI67", "4288",
  "Cell-cycle", "CDKN2A", "1029",
  "Cell-cycle", "CDKN2B", "1030",
  "Cell-cycle", "CDK4", "1019",
  "Cell-cycle", "CDK6", "1021",
  "Cell-cycle", "RB1", "5925",
  "Cell-cycle", "PLK1", "5347",
  "RTK/angiogenesis", "EGFR", "1956",
  "RTK/angiogenesis", "PDGFRA", "5156",
  "RTK/angiogenesis", "KDR", "3791",
  "RTK/angiogenesis", "VEGFA", "7422",
  "Antigen presentation", "B2M", "567",
  "Antigen presentation", "HLA-A", "3105",
  "Antigen presentation", "HLA-B", "3106",
  "Antigen presentation", "HLA-C", "3107",
  "Antigen presentation", "TAP1", "6890",
  "Antigen presentation", "TAP2", "6891",
  "JAK/STAT immune", "JAK1", "3716",
  "JAK/STAT immune", "JAK2", "3717",
  "JAK/STAT immune", "STAT1", "6772",
  "JAK/STAT immune", "STAT3", "6774",
  "Checkpoint context", "CD274", "29126",
  "Checkpoint context", "PDCD1", "5133",
  "Checkpoint context", "HAVCR2", "84868",
  "Checkpoint context", "CD47", "961"
)

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 5 || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(tibble(rho = NA_real_, p_value = NA_real_, n = sum(keep)))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  tibble(rho = unname(test$estimate), p_value = test$p.value, n = sum(keep))
}

message("Running GSE230168 microenvironment and biomarker context analyses")
expr <- read_processed_matrix("data/processed/GSE230168_log2cpm_entrez.csv")
meta <- read_csv("data/processed/GSE230168_metadata_curated.csv", show_col_types = FALSE)
axis_scores <- read_csv("results/GSE230168_axis_scores.csv", show_col_types = FALSE)

common <- intersect(colnames(expr), meta$sample)
expr <- expr[, common, drop = FALSE]
meta <- meta %>% filter(sample %in% common)

micro_scores <- score_gene_sets(expr, microenvironment_sets()) %>%
  left_join(meta, by = "sample")
micro_stats <- compare_scores(
  micro_scores %>% select(sample, starts_with("MICRO_")),
  meta,
  "type",
  "chordoma_I",
  "chordoma_C"
) %>%
  mutate(
    contrast = "chordoma_I_vs_chordoma_C",
    microenvironment_label = micro_labels[pathway]
  ) %>%
  relocate(contrast, microenvironment_label, .after = pathway)

biomarker_expression <- bind_rows(lapply(seq_len(nrow(biomarker_panel)), function(i) {
  gid <- biomarker_panel$entrez_id[i]
  values <- if (gid %in% rownames(expr)) expr[gid, meta$sample] else rep(NA_real_, nrow(meta))
  tibble(
    sample = meta$sample,
    gene_symbol = biomarker_panel$gene_symbol[i],
    entrez_id = gid,
    axis = biomarker_panel$axis[i],
    expression = as.numeric(values)
  )
})) %>%
  left_join(meta, by = "sample")

biomarker_summary <- biomarker_expression %>%
  group_by(axis, gene_symbol, entrez_id, type) %>%
  summarise(
    n = sum(!is.na(expression)),
    mean_log2cpm = mean(expression, na.rm = TRUE),
    median_log2cpm = median(expression, na.rm = TRUE),
    .groups = "drop"
  )

biomarker_group_tests <- biomarker_expression %>%
  filter(type %in% c("chordoma_I", "chordoma_C")) %>%
  group_by(axis, gene_symbol, entrez_id) %>%
  summarise(
    n_chordoma_I = sum(type == "chordoma_I" & !is.na(expression)),
    n_chordoma_C = sum(type == "chordoma_C" & !is.na(expression)),
    mean_chordoma_I = mean(expression[type == "chordoma_I"], na.rm = TRUE),
    mean_chordoma_C = mean(expression[type == "chordoma_C"], na.rm = TRUE),
    delta_chordoma_I_minus_C = mean_chordoma_I - mean_chordoma_C,
    p_value = if (n_chordoma_I >= 2 && n_chordoma_C >= 2) {
      suppressWarnings(wilcox.test(expression[type == "chordoma_I"], expression[type == "chordoma_C"], exact = FALSE)$p.value)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

axis_for_cor <- axis_scores %>%
  select(sample, immune_interferon, cell_cycle, chordoma_marker) %>%
  left_join(micro_scores %>% select(sample, starts_with("MICRO_")), by = "sample")
axis_cols <- setdiff(colnames(axis_for_cor), "sample")
cor_data <- biomarker_expression %>%
  select(sample, gene_symbol, entrez_id, axis, expression, type) %>%
  filter(type %in% c("chordoma_I", "chordoma_C")) %>%
  left_join(axis_for_cor, by = "sample")

biomarker_correlations <- bind_rows(lapply(seq_len(nrow(biomarker_panel)), function(i) {
  sub <- cor_data %>% filter(entrez_id == biomarker_panel$entrez_id[i])
  bind_rows(lapply(axis_cols, function(axis_col) {
    safe_cor(sub$expression, sub[[axis_col]]) %>%
      mutate(
        biomarker_axis = biomarker_panel$axis[i],
        gene_symbol = biomarker_panel$gene_symbol[i],
        entrez_id = biomarker_panel$entrez_id[i],
        correlated_feature = axis_col
      )
  }))
})) %>%
  mutate(
    feature_label = recode(
      correlated_feature,
      immune_interferon = "Immune/interferon axis",
      cell_cycle = "Cell-cycle axis",
      chordoma_marker = "Chordoma marker axis",
      !!!as.list(micro_labels)
    ),
    fdr = p.adjust(p_value, method = "BH")
  ) %>%
  select(biomarker_axis, gene_symbol, entrez_id, correlated_feature, feature_label, n, rho, p_value, fdr)

write_csv(micro_scores, "results/GSE230168_microenvironment_scores.csv")
write_csv(micro_stats, "results/GSE230168_microenvironment_score_tests.csv")
write_csv(biomarker_panel, "results/GSE230168_biomarker_panel.csv")
write_csv(biomarker_expression, "results/GSE230168_biomarker_expression_long.csv")
write_csv(biomarker_summary, "results/GSE230168_biomarker_expression_summary.csv")
write_csv(biomarker_group_tests, "results/GSE230168_biomarker_group_tests.csv")
write_csv(biomarker_correlations, "results/GSE230168_biomarker_axis_correlations.csv")

palette <- group_palette()
sample_order <- micro_scores %>%
  arrange(type, desc(MICRO_T_CELL_INFLAMMATION), sample) %>%
  pull(sample)

micro_long <- micro_scores %>%
  select(sample, type, starts_with("MICRO_")) %>%
  pivot_longer(cols = starts_with("MICRO_"), names_to = "microenvironment", values_to = "score") %>%
  group_by(microenvironment) %>%
  mutate(score_z = as.numeric(scale(score))) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    label = factor(micro_labels[microenvironment], levels = rev(micro_labels))
  )

p_micro <- ggplot(micro_long, aes(x = sample, y = label, fill = score_z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#1E5FA3", mid = "#F1E7DE", high = "#A6492C", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "z-score") +
  labs(title = "a", x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

selected_biomarkers <- c("TBXT", "MKI67", "CDKN2A", "PLK1", "EGFR", "PDGFRA", "KDR", "B2M", "HLA-A", "STAT1", "CD274", "CD47")
biomarker_heat <- biomarker_expression %>%
  filter(gene_symbol %in% selected_biomarkers) %>%
  group_by(gene_symbol) %>%
  mutate(expr_z = as.numeric(scale(expression))) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    gene_symbol = factor(gene_symbol, levels = rev(selected_biomarkers))
  )

p_bio <- ggplot(biomarker_heat, aes(x = sample, y = gene_symbol, fill = expr_z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#1E5FA3", mid = "#F1E7DE", high = "#A6492C", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "z-score") +
  labs(title = "b", x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

top_cor <- biomarker_correlations %>%
  filter(!is.na(rho)) %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = 18) %>%
  mutate(pair = paste(gene_symbol, feature_label, sep = " vs "))

p_cor <- ggplot(top_cor, aes(x = rho, y = reorder(pair, rho), color = fdr < 0.10)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_point(size = 2.3) +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey45"), name = "FDR < 0.10") +
  labs(title = "c", x = "Spearman rho", y = NULL) +
  theme_nature(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

fig <- plot_grid(
  plot_grid(p_micro, p_bio, nrow = 1),
  p_cor,
  ncol = 1,
  rel_heights = c(1, 0.9)
)
ggsave("figures/Supplementary_Figure_GSE230168_biomarkers_microenvironment.png", fig, width = 13, height = 7, dpi = 300)
ggsave("figures/Supplementary_Figure_GSE230168_biomarkers_microenvironment.pdf", fig, width = 13, height = 7)

writeLines(c(
  "# GSE230168 biomarker and microenvironment context",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Processing notes",
  "- Microenvironment values are marker-set scores from bulk RNA-seq, not cell-type deconvolution.",
  "- Biomarker correlations are Spearman correlations in chordoma tumors only.",
  "",
  "## Main outputs",
  "- results/GSE230168_microenvironment_scores.csv",
  "- results/GSE230168_biomarker_axis_correlations.csv",
  "- figures/Supplementary_Figure_GSE230168_biomarkers_microenvironment.png"
), "results/GSE230168_microenvironment_biomarker_summary.md")

message("Done: GSE230168 microenvironment and biomarker context")
