#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Build the primary GSE230168 patient-tumor/nucleus-pulposus landscape.
#
# Main inputs:
#   - data/raw/GSE230168_count_table.csv
#   - data/raw/GSE230168_spl_clust.csv
#   - data/raw/NCBI_gene2ensembl.gz
#   - data/raw/h.all.v2025.1.Hs.entrez.gmt
#
# Main outputs:
#   - data/processed/GSE230168_log2cpm_entrez.csv
#   - results/GSE230168_pathway_scores.csv
#   - results/GSE230168_axis_scores.csv
#   - results/GSE230168_pathway_score_tests.csv
#   - figures/Intermediate_GSE230168_pathway_landscape.*
#
# Statistical approach:
#   Ensembl-to-Entrez harmonization, low-expression filtering, log2(CPM+1),
#   PCA on top-variable genes, mean z-score pathway/signature scoring, and
#   Wilcoxon rank-sum tests with Benjamini-Hochberg FDR for group contrasts.
#
# Interpretation notes:
#   GSE230168 is the anchor cohort. Results define biological axes and
#   signature-level contrasts; they are not used to train a clinical predictor.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

count_path <- "data/raw/GSE230168_count_table.csv"
meta_path <- "data/raw/GSE230168_spl_clust.csv"

message("Reading GSE230168 RNA-seq counts and metadata")
count_df <- read_csv(count_path, show_col_types = FALSE) %>%
  rename(ensembl_gene_id = gene)

meta <- read_csv(meta_path, show_col_types = FALSE) %>%
  select(sample, sex, age, CCJ, days, alive, source, mod1, type) %>%
  mutate(
    type = if_else(is.na(type) & source == "nucleus_pulposus", "nucleus_pulposus", type),
    type = factor(type, levels = c("nucleus_pulposus", "chordoma_I", "chordoma_C")),
    source = factor(source)
  )

sample_cols <- intersect(colnames(count_df)[-1], meta$sample)
missing_from_counts <- setdiff(meta$sample, colnames(count_df)[-1])
if (length(missing_from_counts) > 0) {
  message("Metadata samples not present in count matrix: ", paste(missing_from_counts, collapse = ", "))
}

meta <- meta %>%
  filter(sample %in% sample_cols) %>%
  arrange(match(sample, sample_cols))

gene_map <- read_human_gene2ensembl()
counts_mapped <- map_ensembl_matrix_to_entrez(count_df, "ensembl_gene_id", meta$sample, gene_map)
count_matrix <- matrix_from_gene_table(counts_mapped, "entrez_id", meta$sample)
log_cpm <- log2_cpm(count_matrix, min_cpm = 1, min_samples = 4)

message("Retained ", nrow(log_cpm), " expressed Entrez genes")

gene_sets <- selected_gene_sets()
scores <- score_gene_sets(log_cpm, gene_sets)
axis_scores <- axis_scores_from_pathways(scores)

pathway_scores <- scores %>%
  left_join(meta, by = "sample")

pathway_stats <- bind_rows(
  compare_scores(scores, meta, "type", "chordoma_I", "chordoma_C") %>%
    mutate(contrast = "chordoma_I_vs_chordoma_C"),
  compare_scores(scores, meta, "type", "chordoma_I", "nucleus_pulposus") %>%
    mutate(contrast = "chordoma_I_vs_nucleus_pulposus"),
  compare_scores(scores, meta, "type", "chordoma_C", "nucleus_pulposus") %>%
    mutate(contrast = "chordoma_C_vs_nucleus_pulposus")
) %>%
  relocate(contrast)

axis_scores_annotated <- axis_scores %>%
  left_join(meta, by = "sample")

write_csv(meta, "data/processed/GSE230168_metadata_curated.csv")
write_csv(as.data.frame(log_cpm) %>% rownames_to_column("entrez_id"), "data/processed/GSE230168_log2cpm_entrez.csv")
write_csv(pathway_scores, "results/GSE230168_pathway_scores.csv")
write_csv(gene_set_overlap(log_cpm, gene_sets), "results/GSE230168_pathway_gene_overlap.csv")
write_csv(axis_scores_annotated, "results/GSE230168_axis_scores.csv")
write_csv(pathway_stats, "results/GSE230168_pathway_score_tests.csv")

message("Running PCA")
gene_vars <- apply(log_cpm, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(3000, length(gene_vars)))]
pca <- prcomp(t(log_cpm[top_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_var <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(meta, by = "sample")
write_csv(pca_df, "results/GSE230168_pca_top_variable_genes.csv")

heatmap_sets <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "CUSTOM_ANTIGEN_PRESENTATION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "CUSTOM_CDK46_AXIS",
  "CUSTOM_PLK1_AXIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS",
  "CUSTOM_DRUGGABLE_RTKS",
  "CUSTOM_CHORDOMA_MARKERS"
)

labels <- pathway_labels()
sample_order <- meta$sample[order(meta$type, meta$sample)]
palette <- group_palette()

score_long <- pathway_scores %>%
  select(sample, type, all_of(heatmap_sets)) %>%
  pivot_longer(cols = all_of(heatmap_sets), names_to = "pathway", values_to = "score") %>%
  group_by(pathway) %>%
  mutate(score_z = as.numeric(scale(score))) %>%
  ungroup() %>%
  mutate(
    pathway_label = factor(labels[pathway], levels = rev(labels[heatmap_sets])),
    sample = factor(sample, levels = sample_order)
  )

group_counts <- meta %>% count(type)
sample_annotation <- meta %>%
  mutate(sample = factor(sample, levels = sample_order), annotation = "Group")

p_design <- ggplot(group_counts, aes(x = type, y = n, fill = type)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_text(aes(label = n), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(title = "A. GSE230168 study groups", x = NULL, y = "Samples") +
  coord_cartesian(ylim = c(0, max(group_counts$n) + 3)) +
  theme_nature(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(title = "B. RNA-seq PCA", x = paste0("PC1 (", pca_var[1], "%)"), y = paste0("PC2 (", pca_var[2], "%)"), color = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

p_axes <- ggplot(axis_scores_annotated, aes(x = immune_interferon, y = cell_cycle, color = type)) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(title = "C. Translational axes", x = "Immune/interferon score", y = "Cell-cycle score", color = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

p_heat_annot <- ggplot(sample_annotation, aes(x = sample, y = annotation, fill = type)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = palette, drop = FALSE, name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "none", plot.margin = margin(0, 6, 0, 6))

p_heat <- ggplot(score_long, aes(x = sample, y = pathway_label, fill = score_z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#1E5FA3", mid = "#F1E7DE", high = "#A6492C", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "z-score") +
  labs(title = "D. Pathway/signature scores", x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), legend.position = "right", plot.margin = margin(0, 6, 6, 6))

heat_block <- plot_grid(p_heat_annot, p_heat, ncol = 1, rel_heights = c(0.04, 1), align = "v", axis = "lr")
fig <- plot_grid(plot_grid(p_design, p_pca, p_axes, nrow = 1, rel_widths = c(0.9, 1.15, 1.15)), heat_block, ncol = 1, rel_heights = c(1, 1.25))

ggsave("figures/Intermediate_GSE230168_pathway_landscape.png", fig, width = 13, height = 8.5, dpi = 300)
ggsave("figures/Intermediate_GSE230168_pathway_landscape.pdf", fig, width = 13, height = 8.5)

candidate_genes <- tibble(
  axis = c("TBXT/chordoma identity", "JAK/STAT immune", "JAK/STAT immune", "Cell-cycle/CDK4-6", "Cell-cycle/CDK4-6", "Cell-cycle/PLK1", "RTK", "RTK", "RTK/angiogenesis", "Antigen presentation"),
  gene_symbol = c("TBXT", "JAK1", "STAT3", "CDK4", "CDKN2A", "PLK1", "EGFR", "PDGFRA", "KDR", "B2M"),
  entrez_id = c("6862", "3716", "6774", "1019", "1029", "5347", "1956", "5156", "3791", "567")
)

gene_summary <- lapply(seq_len(nrow(candidate_genes)), function(i) {
  gid <- candidate_genes$entrez_id[i]
  values <- if (gid %in% rownames(log_cpm)) log_cpm[gid, meta$sample] else rep(NA_real_, nrow(meta))
  tibble(sample = meta$sample, value = as.numeric(values), type = meta$type) %>%
    group_by(type) %>%
    summarise(mean_log2cpm = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(gene_symbol = candidate_genes$gene_symbol[i], entrez_id = gid)
}) %>%
  bind_rows() %>%
  left_join(candidate_genes, by = c("gene_symbol", "entrez_id")) %>%
  select(axis, gene_symbol, entrez_id, type, mean_log2cpm)

write_csv(gene_summary, "results/GSE230168_candidate_gene_expression.csv")

summary_tbl <- pathway_scores %>%
  select(sample, type, all_of(heatmap_sets)) %>%
  pivot_longer(cols = all_of(heatmap_sets), names_to = "pathway", values_to = "score") %>%
  group_by(type, pathway) %>%
  summarise(n = n(), mean_score = mean(score, na.rm = TRUE), sd_score = sd(score, na.rm = TRUE), .groups = "drop")
write_csv(summary_tbl, "results/GSE230168_pathway_scores_by_group.csv")

writeLines(c(
  "# GSE230168 first-pass pathway analysis",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Processing notes",
  paste0("- Count matrix input rows: ", nrow(count_df), "."),
  paste0("- Final groups: ", paste(group_counts$type, group_counts$n, sep = "=", collapse = "; "), "."),
  paste0("- Expressed Entrez genes retained after CPM filtering: ", nrow(log_cpm), "."),
  "- Pathway/signature scores were calculated as mean sample-wise z-scores.",
  "- Wilcoxon tests with Benjamini-Hochberg FDR were used for group comparisons.",
  "",
  "## Main outputs",
  "- figures/Intermediate_GSE230168_pathway_landscape.png",
  "- results/GSE230168_pathway_score_tests.csv",
  "- results/GSE230168_candidate_gene_expression.csv"
), "results/GSE230168_analysis_summary.md")

write_session_info()
message("Done: GSE230168 Figure 1 and tables")
