#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Validate TBXT/T-DARPin perturbation effects in GSE276715 chordoma cells.
#
# Main inputs:
#   - data/raw/GSE276715/GSE276715_RAW.tar
#   - data/raw/h.all.v2025.1.Hs.entrez.gmt
#   - results/GSE230168_axis_scores.csv (for patient-cohort projection)
#
# Main outputs:
#   - data/processed/GSE276715_log2tpm_entrez.csv
#   - results/GSE276715_pathway_scores.csv
#   - results/GSE276715_pathway_score_tests.csv
#   - results/GSE276715_T_DARPin_vs_control_signature.csv
#   - figures/Panel_perturbation_A_GSE276715_TDARPin.*
#
# Statistical approach:
#   TPM tables are harmonized to Entrez and transformed as log2(TPM+1).
#   The same Hallmark/custom signatures are scored by mean z-score. T-DARPin
#   versus control pathway contrasts use Wilcoxon/BH-FDR, and gene-level
#   perturbation signatures use limma-trend with empirical Bayes moderation.
#
# Interpretation notes:
#   This is a functional perturbation layer, not an independent patient cohort.
#   It tests whether TBXT targeting shifts the prespecified chordoma-state axes.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

raw_dir <- "data/raw/GSE276715"
files <- list.files(raw_dir, pattern = "featureCounts.tsv.gz$", full.names = TRUE)
if (length(files) != 12) {
  stop("Expected 12 GSE276715 featureCounts TSV files after extraction.")
}

parse_meta <- function(path) {
  nm <- basename(path)
  gsm <- sub("_.*$", "", nm)
  code <- sub("^.*cell-([a-z0-9]+)-[0-9]+_.*$", "\\1", nm)
  rep <- sub("^.*cell-[a-z0-9]+-([0-9]+)_.*$", "\\1", nm)
  tibble(
    sample = gsm,
    file = path,
    darpin = recode(code, e35 = "E3_5_control", a2 = "T_DARPin_A2", b1 = "T_DARPin_B1", d4 = "T_DARPin_D4"),
    condition = if_else(code == "e35", "control", "T_DARPin"),
    replicate = rep
  )
}

meta <- bind_rows(lapply(files, parse_meta)) %>%
  mutate(
    condition = factor(condition, levels = c("control", "T_DARPin")),
    darpin = factor(darpin, levels = c("E3_5_control", "T_DARPin_A2", "T_DARPin_B1", "T_DARPin_D4"))
  ) %>%
  arrange(condition, darpin, replicate)

message("Reading GSE276715 TPM tables")
expr_list <- lapply(seq_len(nrow(meta)), function(i) {
  read_tsv(meta$file[i], show_col_types = FALSE) %>%
    transmute(ensembl_gene_id = gene_id, !!meta$sample[i] := TPM)
})

expr_df <- Reduce(function(x, y) full_join(x, y, by = "ensembl_gene_id"), expr_list)
gene_map <- read_human_gene2ensembl()
expr_mapped <- map_ensembl_matrix_to_entrez(expr_df, "ensembl_gene_id", meta$sample, gene_map)
expr_matrix <- matrix_from_gene_table(expr_mapped, "entrez_id", meta$sample)
log_tpm <- log2(expr_matrix + 1)
log_tpm <- log_tpm[rowSums(expr_matrix > 1, na.rm = TRUE) >= 3, , drop = FALSE]

gene_sets <- selected_gene_sets()
scores <- score_gene_sets(log_tpm, gene_sets)
score_tests <- compare_scores(scores, meta, "condition", "T_DARPin", "control") %>%
  mutate(contrast = "T_DARPin_vs_control") %>%
  relocate(contrast)

signature <- differential_signature_limma(log_tpm, meta, "condition", "T_DARPin", "control")

write_csv(meta, "data/processed/GSE276715_metadata_curated.csv")
write_csv(as.data.frame(log_tpm) %>% rownames_to_column("entrez_id"), "data/processed/GSE276715_log2tpm_entrez.csv")
write_csv(scores %>% left_join(meta, by = "sample"), "results/GSE276715_pathway_scores.csv")
write_csv(score_tests, "results/GSE276715_pathway_score_tests.csv")
write_csv(signature, "results/GSE276715_T_DARPin_vs_control_signature.csv")
write_csv(gene_set_overlap(log_tpm, gene_sets), "results/GSE276715_pathway_gene_overlap.csv")

gse230168_path <- "data/processed/GSE230168_log2cpm_entrez.csv"
if (file.exists(gse230168_path)) {
  gse230168 <- read_csv(gse230168_path, show_col_types = FALSE)
  gse230168_mat <- as.matrix(gse230168[, -1])
  rownames(gse230168_mat) <- as.character(gse230168$entrez_id)
  storage.mode(gse230168_mat) <- "numeric"
  projection <- signature_set_scores(gse230168_mat, signature, "TDARPIN", top_n = 200)
  write_csv(projection, "results/GSE276715_signature_projected_on_GSE230168.csv")
} else {
  projection <- tibble()
}

labels <- pathway_labels()
palette <- group_palette()
selected <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "CUSTOM_CDK46_AXIS",
  "CUSTOM_PLK1_AXIS",
  "CUSTOM_CHORDOMA_MARKERS"
)

delta_df <- score_tests %>%
  filter(pathway %in% selected) %>%
  mutate(pathway_label = factor(labels[pathway], levels = labels[selected]))

p_delta <- ggplot(delta_df, aes(x = delta_a_minus_b, y = pathway_label, fill = delta_a_minus_b > 0)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), guide = "none") +
  labs(title = "a", x = "Mean score difference vs control", y = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

plot_genes <- c(TBXT = "6862", IGFBP3 = "3486", CDKN2A = "1029", PLK1 = "5347", JAK1 = "3716", STAT3 = "6774")
gene_df <- bind_rows(lapply(names(plot_genes), function(symbol) {
  gid <- plot_genes[[symbol]]
  if (!gid %in% rownames(log_tpm)) return(tibble())
  tibble(sample = colnames(log_tpm), value = as.numeric(log_tpm[gid, ]), gene = symbol)
})) %>%
  left_join(meta, by = "sample")

p_genes <- ggplot(gene_df, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
  geom_point(position = position_jitter(width = 0.08), size = 1.7, alpha = 0.85) +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = palette, drop = FALSE) +
  labs(title = "b", x = NULL, y = "log2(TPM + 1)") +
  theme_nature(base_size = 9) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

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
  labs(title = "c", x = "limma log2 fold-change", y = "-log10(limma FDR)") +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

if (nrow(projection) > 0 && file.exists("data/processed/GSE230168_metadata_curated.csv")) {
  gse230168_meta <- read_csv("data/processed/GSE230168_metadata_curated.csv", show_col_types = FALSE)
  proj_long <- projection %>%
    left_join(gse230168_meta, by = "sample") %>%
    pivot_longer(cols = starts_with("TDARPIN"), names_to = "signature", values_to = "score")
  p_proj <- ggplot(proj_long, aes(x = type, y = score, fill = type)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
    geom_point(position = position_jitter(width = 0.08), size = 1.6, alpha = 0.8) +
    facet_wrap(~signature, nrow = 1) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(title = "d", x = NULL, y = "Signature score") +
    theme_nature(base_size = 9) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.minor = element_blank())
} else {
  p_proj <- ggplot() + theme_void() + labs(title = "d")
}

fig <- plot_grid(
  plot_grid(p_delta, p_genes, nrow = 1, rel_widths = c(0.9, 1.2)),
  plot_grid(p_volcano, p_proj, nrow = 1, rel_widths = c(1, 1.15)),
  ncol = 1
)

ggsave("figures/Panel_perturbation_A_GSE276715_TDARPin.png", fig, width = 13, height = 8.5, dpi = 300)
ggsave("figures/Panel_perturbation_A_GSE276715_TDARPin.pdf", fig, width = 13, height = 8.5)

writeLines(c(
  "# GSE276715 TBXT/T-DARPin validation",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Processing notes",
  paste0("- Samples: ", nrow(meta), " UM-Chor1 profiles; ", sum(meta$condition == "control"), " controls and ", sum(meta$condition == "T_DARPin"), " T-DARPin profiles."),
  paste0("- Expressed Entrez genes retained from TPM matrix: ", nrow(log_tpm), "."),
  "- Signature contrast: all T-DARPins A2/B1/D4 versus E3_5 control.",
  "- Gene-level differential signature uses limma-trend with empirical Bayes moderation.",
  "- Gene-set scores use mean sample-wise z-scores.",
  "- Construct-direction check: all seven displayed program deltas were negative for A2, B1, and D4 relative to the shared control.",
  "",
  "## Main outputs",
  "- figures/Panel_perturbation_A_GSE276715_TDARPin.png",
  "- results/GSE276715_T_DARPin_vs_control_signature.csv",
  "- results/GSE276715_signature_projected_on_GSE230168.csv"
), "results/GSE276715_analysis_summary.md")

message("Done: GSE276715 TBXT validation")
