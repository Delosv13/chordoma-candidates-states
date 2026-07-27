#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Validate CDK4/6-palbociclib response biology in GSE216418 PDX models.
#
# Main inputs:
#   - data/raw/GSE216418/GSE216418_RAW.tar
#   - data/raw/NCBI_gene2ensembl.gz
#   - data/raw/h.all.v2025.1.Hs.entrez.gmt
#
# Main outputs:
#   - data/processed/GSE216418_log2cpm_entrez.csv
#   - results/GSE216418_pathway_scores.csv
#   - results/GSE216418_axis_scores.csv
#   - results/GSE216418_palbociclib_vs_untreated_signature.csv
#   - figures/Panel_perturbation_B_GSE216418_palbociclib.*
#
# Statistical approach:
#   Count tables are mapped to Entrez, filtered, and transformed to
#   log2(CPM+1). Pathway/axis contrasts use Wilcoxon/BH-FDR. Gene-level
#   palbociclib signatures use limma-trend/eBayes with PDX model blocking.
#
# Interpretation notes:
#   Blocking by PDX model reduces confounding between model background and
#   treatment effect. This layer supports mechanism, not clinical response
#   prediction.
#
# Panel lettering (2026-07-26):
#   Sub-panel titles are lettered e-h, not a-d, because this figure is the
#   lower half of the merged main-text Figure 3. Script 03 supplies panels
#   a-d (GSE276715 / T-DARPin) and this script supplies panels e-h
#   (GSE216418 / palbociclib), so the assembled figure carries one
#   uninterrupted a-h sequence and no longer needs dataset banners.
#   See scripts/extras/merge_perturbation_figure.R.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

raw_dir <- "data/raw/GSE216418"
files <- list.files(raw_dir, pattern = "read_counts.txt.gz$", full.names = TRUE)
if (length(files) != 12) {
  stop("Expected 12 GSE216418 read count files after extraction.")
}

meta_lookup <- tribble(
  ~sample, ~model, ~treatment, ~replicate,
  "GSM6672200", "CD3", "untreated", "1",
  "GSM6672201", "CD3", "untreated", "2",
  "GSM6672202", "CD3", "untreated", "3",
  "GSM6672203", "CD7", "untreated", "1",
  "GSM6672204", "CD7", "untreated", "2",
  "GSM6672205", "CD7", "untreated", "3",
  "GSM6672206", "CD3", "palbociclib", "1",
  "GSM6672207", "CD3", "palbociclib", "2",
  "GSM6672208", "CD3", "palbociclib", "3",
  "GSM6672209", "CD7", "palbociclib", "1",
  "GSM6672210", "CD7", "palbociclib", "2",
  "GSM6672211", "CD7", "palbociclib", "3"
)

file_meta <- tibble(
  sample = sub("_.*$", "", basename(files)),
  file = files
) %>%
  left_join(meta_lookup, by = "sample") %>%
  mutate(treatment = factor(treatment, levels = c("untreated", "palbociclib")))

if (any(is.na(file_meta$treatment))) {
  stop("Failed to map one or more GSE216418 files to metadata.")
}

file_meta <- file_meta %>% arrange(model, treatment, replicate)

message("Reading GSE216418 read-count tables")
count_list <- lapply(seq_len(nrow(file_meta)), function(i) {
  read_tsv(file_meta$file[i], col_names = c("ensembl_gene_id", file_meta$sample[i]), show_col_types = FALSE)
})

count_df <- Reduce(function(x, y) full_join(x, y, by = "ensembl_gene_id"), count_list)
gene_map <- read_human_gene2ensembl()
counts_mapped <- map_ensembl_matrix_to_entrez(count_df, "ensembl_gene_id", file_meta$sample, gene_map)
count_matrix <- matrix_from_gene_table(counts_mapped, "entrez_id", file_meta$sample)
log_cpm <- log2_cpm(count_matrix, min_cpm = 1, min_samples = 3)

gene_sets <- selected_gene_sets()
scores <- score_gene_sets(log_cpm, gene_sets)
score_tests <- compare_scores(scores, file_meta, "treatment", "palbociclib", "untreated") %>%
  mutate(contrast = "palbociclib_vs_untreated") %>%
  relocate(contrast)
signature <- differential_signature_limma(log_cpm, file_meta, "treatment", "palbociclib", "untreated", block_col = "model")

axis_scores <- axis_scores_from_pathways(scores) %>%
  left_join(file_meta, by = "sample")

write_csv(file_meta, "data/processed/GSE216418_metadata_curated.csv")
write_csv(as.data.frame(log_cpm) %>% rownames_to_column("entrez_id"), "data/processed/GSE216418_log2cpm_entrez.csv")
write_csv(scores %>% left_join(file_meta, by = "sample"), "results/GSE216418_pathway_scores.csv")
write_csv(score_tests, "results/GSE216418_pathway_score_tests.csv")
write_csv(axis_scores, "results/GSE216418_axis_scores.csv")
write_csv(signature, "results/GSE216418_palbociclib_vs_untreated_signature.csv")
write_csv(gene_set_overlap(log_cpm, gene_sets), "results/GSE216418_pathway_gene_overlap.csv")

labels <- pathway_labels()
palette <- group_palette()
selected <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "CUSTOM_CDK46_AXIS",
  "CUSTOM_PLK1_AXIS",
  "HALLMARK_MTORC1_SIGNALING",
  "CUSTOM_CHORDOMA_MARKERS"
)

delta_df <- score_tests %>%
  filter(pathway %in% selected) %>%
  mutate(pathway_label = factor(labels[pathway], levels = labels[selected]))

p_delta <- ggplot(delta_df, aes(x = delta_a_minus_b, y = pathway_label, fill = delta_a_minus_b > 0)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), guide = "none") +
  labs(title = "e", x = "Mean score difference vs untreated", y = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

axis_long <- axis_scores %>%
  select(sample, model, treatment, immune_interferon, cell_cycle, chordoma_marker) %>%
  pivot_longer(cols = c(immune_interferon, cell_cycle, chordoma_marker), names_to = "axis", values_to = "score") %>%
  mutate(axis = recode(axis, immune_interferon = "Immune/interferon", cell_cycle = "Cell cycle", chordoma_marker = "Chordoma marker"))

p_axis <- ggplot(axis_long, aes(x = treatment, y = score, fill = treatment)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
  geom_point(aes(shape = model), position = position_jitter(width = 0.08), size = 1.8, alpha = 0.9) +
  facet_wrap(~axis, nrow = 1) +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(title = "f", x = NULL, y = "Score", shape = "PDX") +
  theme_nature(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())

plot_genes <- c(CDK4 = "1019", CDK6 = "1021", RB1 = "5925", E2F1 = "1869", PLK1 = "5347", CDKN2A = "1029")
gene_df <- bind_rows(lapply(names(plot_genes), function(symbol) {
  gid <- plot_genes[[symbol]]
  if (!gid %in% rownames(log_cpm)) return(tibble())
  tibble(sample = colnames(log_cpm), value = as.numeric(log_cpm[gid, ]), gene = symbol)
})) %>%
  left_join(file_meta, by = "sample")

p_genes <- ggplot(gene_df, aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
  geom_point(aes(shape = model), position = position_jitter(width = 0.08), size = 1.6, alpha = 0.85) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(title = "g", x = NULL, y = "log2(CPM + 1)", shape = "PDX") +
  theme_nature(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())

volcano_df <- signature %>%
  mutate(
    neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
    label = if_else(entrez_id %in% plot_genes, names(plot_genes)[match(entrez_id, plot_genes)], NA_character_)
  )

p_volcano <- ggplot(volcano_df, aes(x = log2fc, y = neg_log10_fdr)) +
  geom_point(color = "grey60", size = 0.8, alpha = 0.5) +
  geom_point(data = filter(volcano_df, !is.na(label)), color = "#D55E00", size = 2) +
  ggrepel::geom_text_repel(aes(label = label), max.overlaps = 20, size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  labs(title = "h", x = "limma log2 fold-change", y = "-log10(limma FDR)") +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

fig <- plot_grid(
  plot_grid(p_delta, p_axis, nrow = 1, rel_widths = c(0.9, 1.25)),
  plot_grid(p_genes, p_volcano, nrow = 1, rel_widths = c(1.2, 1)),
  ncol = 1
)

ggsave("figures/Panel_perturbation_B_GSE216418_palbociclib.png", fig, width = 13, height = 8.5, dpi = 300)
ggsave("figures/Panel_perturbation_B_GSE216418_palbociclib.pdf", fig, width = 13, height = 8.5)

writeLines(c(
  "# GSE216418 palbociclib validation",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Processing notes",
  paste0("- Samples: ", nrow(file_meta), " PDX RNA-seq profiles across CD3 and CD7."),
  paste0("- Expressed Entrez genes retained after CPM filtering: ", nrow(log_cpm), "."),
  "- Signature contrast: palbociclib versus untreated across both PDX models.",
  "- Gene-level differential signature uses limma-trend with empirical Bayes moderation and PDX model as a blocking factor.",
  "- Gene-set scores use mean sample-wise z-scores.",
  "",
  "## Main outputs",
  "- figures/Panel_perturbation_B_GSE216418_palbociclib.png",
  "- results/GSE216418_palbociclib_vs_untreated_signature.csv",
  "- results/GSE216418_pathway_score_tests.csv"
), "results/GSE216418_analysis_summary.md")

message("Done: GSE216418 palbociclib validation")
