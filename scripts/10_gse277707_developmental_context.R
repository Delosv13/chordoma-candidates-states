#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Add developmental/TBXT-rs2305089 context from GSE277707 differentiation.
#
# Main inputs:
#   - data/raw/GSE277707/GSE277707_norm_matrix_48and72hrs.csv.gz
#   - data/raw/h.all.v2025.1.Hs.entrez.gmt
#
# Main outputs:
#   - data/processed/GSE277707_log2norm_entrez.csv
#   - results/GSE277707_pathway_scores.csv
#   - results/GSE277707_pathway_score_tests_all.csv
#   - results/GSE277707_pathway_score_tests_by_stage.csv
#   - results/GSE277707_heterozygous_vs_wild_type_signature.csv
#   - figures/Supplementary_Figure_GSE277707_rs2305089_notochord_validation.*
#
# Statistical approach:
#   Normalized Entrez-level expression is transformed as log2(value+1).
#   Prespecified signatures are scored by mean z-score and compared by
#   Wilcoxon/BH-FDR overall and within stage. The exploratory gene-level
#   genotype signature uses the legacy Welch t-test helper.
#
# Interpretation notes:
#   This dataset is developmental context, not tumor or treatment validation.
#   Interpret genotype/stage effects as biological plausibility evidence.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

matrix_path <- "data/raw/GSE277707/GSE277707_norm_matrix_48and72hrs.csv.gz"
if (!file.exists(matrix_path)) {
  stop("Missing GSE277707 normalized matrix.")
}

message("Reading GSE277707 normalized matrix")
df <- read_csv(matrix_path, show_col_types = FALSE)
sample_cols <- setdiff(colnames(df), c("ENTREZID", "SYMBOL", "GENENAME"))

meta <- tibble(sample = sample_cols) %>%
  mutate(
    genotype = if_else(startsWith(sample, "463N"), "heterozygous", "wild_type"),
    line = sub("^(4610N|463N)([0-9]+)_.*$", "\\1\\2", sample),
    time_h = as.integer(sub("^.*_([0-9]+)_[0-9]+$", "\\1", sample)),
    replicate = sub("^.*_[0-9]+_([0-9]+)$", "\\1", sample),
    stage = if_else(time_h == 48, "mesoderm_48h", "MSC_72h"),
    genotype = factor(genotype, levels = c("wild_type", "heterozygous")),
    stage = factor(stage, levels = c("mesoderm_48h", "MSC_72h"))
  )

expr_tbl <- df %>%
  transmute(entrez_id = as.character(ENTREZID), across(all_of(sample_cols), as.numeric)) %>%
  filter(!is.na(entrez_id), entrez_id != "NA") %>%
  group_by(entrez_id) %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop")

expr <- matrix_from_gene_table(expr_tbl, "entrez_id", sample_cols)
log_expr <- log2(expr + 1)
log_expr <- log_expr[rowSums(expr > 1, na.rm = TRUE) >= 3, , drop = FALSE]

gene_sets <- selected_gene_sets()
scores <- score_gene_sets(log_expr, gene_sets)
score_tests_all <- compare_scores(scores, meta, "genotype", "heterozygous", "wild_type") %>%
  mutate(contrast = "heterozygous_vs_wild_type_all") %>%
  relocate(contrast)

score_tests_by_stage <- bind_rows(lapply(levels(meta$stage), function(st) {
  meta_st <- meta %>% filter(stage == st)
  compare_scores(scores %>% select(sample, everything()), meta_st, "genotype", "heterozygous", "wild_type") %>%
    mutate(contrast = paste0("heterozygous_vs_wild_type_", st))
})) %>%
  relocate(contrast)

signature_all <- differential_signature(log_expr, meta, "genotype", "heterozygous", "wild_type")

write_csv(meta, "data/processed/GSE277707_metadata_curated.csv")
write_csv(as.data.frame(log_expr) %>% rownames_to_column("entrez_id"), "data/processed/GSE277707_log2norm_entrez.csv")
write_csv(scores %>% left_join(meta, by = "sample"), "results/GSE277707_pathway_scores.csv")
write_csv(score_tests_all, "results/GSE277707_pathway_score_tests_all.csv")
write_csv(score_tests_by_stage, "results/GSE277707_pathway_score_tests_by_stage.csv")
write_csv(signature_all, "results/GSE277707_heterozygous_vs_wild_type_signature.csv")
write_csv(gene_set_overlap(log_expr, gene_sets), "results/GSE277707_pathway_gene_overlap.csv")

labels <- pathway_labels()
palette <- group_palette()
selected <- c(
  "CUSTOM_CHORDOMA_MARKERS",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_E2F_TARGETS"
)

delta_df <- score_tests_by_stage %>%
  filter(pathway %in% selected) %>%
  mutate(
    pathway_label = factor(labels[pathway], levels = labels[selected]),
    stage = sub("^heterozygous_vs_wild_type_", "", contrast)
  )

p_delta <- ggplot(delta_df, aes(x = delta_a_minus_b, y = pathway_label, fill = stage)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, color = "grey25", linewidth = 0.15) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = c(mesoderm_48h = "#009E73", MSC_72h = "#D55E00")) +
  labs(title = "a", x = "Mean score difference vs wild type", y = NULL, fill = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

plot_genes <- c(TBXT = "6862", SOX9 = "6662", WNT5A = "7474", IGFBP3 = "3486", MTOR = "2475", CDK4 = "1019")
gene_df <- bind_rows(lapply(names(plot_genes), function(symbol) {
  gid <- plot_genes[[symbol]]
  if (!gid %in% rownames(log_expr)) return(tibble())
  tibble(sample = colnames(log_expr), value = as.numeric(log_expr[gid, ]), gene = symbol)
})) %>%
  left_join(meta, by = "sample")

p_genes <- ggplot(gene_df, aes(x = genotype, y = value, fill = genotype)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.08), size = 1.5, alpha = 0.85) +
  facet_grid(stage ~ gene, scales = "free_y") +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(title = "b", x = NULL, y = "log2(normalized expression + 1)") +
  theme_nature(base_size = 8.5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())

gene_vars <- apply(log_expr, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(2500, length(gene_vars)))]
pca <- prcomp(t(log_expr[top_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
pca_var <- round(100 * summary(pca)$importance[2, 1:2], 1)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(meta, by = "sample")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = genotype, shape = stage)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(title = "c", x = paste0("PC1 (", pca_var[1], "%)"), y = paste0("PC2 (", pca_var[2], "%)"), color = NULL, shape = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

volcano_df <- signature_all %>%
  mutate(
    neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
    label = if_else(entrez_id %in% plot_genes, names(plot_genes)[match(entrez_id, plot_genes)], NA_character_)
  )

p_volcano <- ggplot(volcano_df, aes(x = log2fc, y = neg_log10_fdr)) +
  geom_point(color = "grey60", size = 0.8, alpha = 0.5) +
  geom_point(data = filter(volcano_df, !is.na(label)), color = "#D55E00", size = 2) +
  ggrepel::geom_text_repel(aes(label = label), max.overlaps = 20, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  labs(title = "d", x = "log2 fold-change", y = "-log10(FDR)") +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

fig <- plot_grid(
  plot_grid(p_delta, p_pca, nrow = 1, rel_widths = c(1, 1)),
  plot_grid(p_genes, p_volcano, nrow = 1, rel_widths = c(1.25, 1)),
  ncol = 1
)

ggsave("figures/Supplementary_Figure_GSE277707_rs2305089_notochord_validation.png", fig, width = 13, height = 8.5, dpi = 300)
ggsave("figures/Supplementary_Figure_GSE277707_rs2305089_notochord_validation.pdf", fig, width = 13, height = 8.5)

writeLines(c(
  "# GSE277707 rs2305089 notochordal validation",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Processing notes",
  paste0("- Samples: ", nrow(meta), " engineered differentiation profiles."),
  paste0("- Expressed Entrez genes retained: ", nrow(log_expr), "."),
  "- Signature contrast: rs2305089 heterozygous versus wild type across 48h and 72h models.",
  "- Stage-specific pathway tests were exported separately.",
  "",
  "## Main outputs",
  "- figures/Supplementary_Figure_GSE277707_rs2305089_notochord_validation.png",
  "- results/GSE277707_pathway_score_tests_by_stage.csv",
  "- results/GSE277707_heterozygous_vs_wild_type_signature.csv"
), "results/GSE277707_analysis_summary.md")

message("Done: GSE277707 rs2305089 validation")
