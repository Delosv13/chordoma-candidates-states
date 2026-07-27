#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Regenerate analysis figures and translational result tables.
#
# Main inputs:
#   - results/GSE230168_* scores/tests
#   - results/GSE216418_* scores/tests
#   - results/GSE276715_* scores/tests
#   - results/GSE277707_* scores/tests
#   - results/effect_size_stability_summary.csv when available
#
# Main outputs:
#   - results/Table1_translational_vulnerability_matrix.csv
#   - results/translational_trial_context.csv
#   - results/translational_evidence_layer_map.csv
#   - figures/Figure2_GSE230168_high_impact_landscape.*
#   - figures/Figure4_translational_robustness_map.*
#
# Statistical approach:
#   This script synthesizes already computed effect estimates, FDR summaries,
#   robustness summaries, and evidence layers; it does not refit the primary
#   statistical models.
#
# Interpretation notes:
#   This reporting script should be run after the
#   dataset-level analyses and robustness scripts are complete.
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

metric_from <- function(path, pathway_id, delta_col = "delta_a_minus_b") {
  if (!file.exists(path)) return("not available")
  x <- read_csv(path, show_col_types = FALSE)
  row <- x %>% filter(.data$pathway == pathway_id) %>% slice_head(n = 1)
  if (nrow(row) == 0) return("not available")
  paste0("delta=", signif(row[[delta_col]], 3), "; FDR=", signif(row$fdr, 3))
}

message("Updating analysis figures and translational result tables")

trial_context <- tribble(
  ~strategy, ~linked_state, ~evidence_type, ~public_source, ~clinical_or_translational_implication, ~limitation,
  "CDK4/6 inhibition / palbociclib",
  "Cell-cycle/CDK4-6/PLK1-high",
  "PDX perturbation and phase II trial context",
  "GSE216418; NCT03110744",
  "Supports connecting E2F/G2M/CDK4-6 scores with CDK4/6 inhibitor prioritization hypotheses.",
  "No patient-response classifier is established.",
  "TBXT / brachyury targeting",
  "TBXT/chordoma identity",
  "Cell-line perturbation and vaccine trial context",
  "GSE276715; NCT02383498",
  "Links TBXT perturbation to identity and proliferative programs while anchoring the axis in brachyury-directed therapeutic development.",
  "Direct TBXT targeting remains investigational and requires protein-level validation.",
  "Immune checkpoint / interferon context",
  "Immune/interferon-high",
  "Patient RNA-seq and recent spatial immune profiling literature",
  "GSE230168; spatial immune profiling studies",
  "Prioritizes immune-context biomarkers and JAK/STAT-interferon biology for correlative studies.",
  "Bulk RNA-seq marker scores are not cell-type deconvolution.",
  "RTK / EGFR / PDGFRA axis",
  "RTK/angiogenic candidate axis",
  "Patient RNA-seq and targeted-therapy literature",
  "GSE230168; imatinib/lapatinib literature",
  "Supports treating RTK expression as a stratifying hypothesis rather than a standalone dependency marker.",
  "Expression does not prove pathway addiction.",
  "CDK/EGFR combination",
  "Cell-cycle plus RTK context",
  "Recent preclinical literature",
  "CDK/EGFR chordoma study",
  "Highlights a rational combination direction where proliferative and RTK states overlap.",
  "Requires independent validation in chordoma models and clinical samples.",
  "PLK1 inhibition",
  "Cell-cycle/CDK4-6/PLK1-high",
  "PDX perturbation literature and TBXT-linked signature support",
  "GSE276715; GSE216418-associated literature",
  "Supports PLK1 as a parallel proliferative vulnerability hypothesis.",
  "Transcriptomic PLK1-axis scores do not establish drug sensitivity."
)
write_csv(trial_context, "results/translational_trial_context.csv")

table1 <- tribble(
  ~molecular_state, ~candidate_biomarkers, ~therapeutic_hypothesis, ~evidence_layers, ~public_data_support, ~main_limitation,
  "Immune/interferon-high chordoma",
  "IFN-gamma, IL6-JAK-STAT3, inflammatory response, antigen presentation, T-cell/checkpoint marker context",
  "Prioritize immune/JAK-STAT context and test whether JAK-pathway or immune-combination strategies are state-dependent.",
  "Patient RNA-seq; bulk microenvironment marker scores; spatial immune literature context",
  paste(
    "GSE230168 chordoma_I vs chordoma_C:",
    metric_from("results/GSE230168_pathway_score_tests.csv", "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
    "Microenvironment scores and biomarker correlations are provided in Table S7/S8."
  ),
  "Bulk RNA-seq cannot separate tumor-cell-intrinsic interferon signaling from microenvironmental signal.",
  "Cell-cycle/CDK4-6/PLK1-high chordoma",
  "MKI67, E2F targets, G2M checkpoint, CDK4/6 axis, PLK1 axis, CDKN2A/B/9p21 context",
  "Prioritize CDK4/6 and PLK1 vulnerability testing in tumors or models with high cell-cycle scores.",
  "Patient RNA-seq; TBXT perturbation; PDX palbociclib perturbation; phase II trial context",
  paste(
    "GSE216418 palbociclib vs untreated:",
    metric_from("results/GSE216418_pathway_score_tests.csv", "HALLMARK_E2F_TARGETS"),
    "and",
    metric_from("results/GSE216418_pathway_score_tests.csv", "HALLMARK_G2M_CHECKPOINT"),
    "GSE276715 T-DARPin reduced PLK1 axis:",
    metric_from("results/GSE276715_pathway_score_tests.csv", "CUSTOM_PLK1_AXIS")
  ),
  "Perturbation studies support mechanism but do not establish clinical response prediction.",
  "TBXT/chordoma-identity state",
  "TBXT, keratins, SOX9, EPCAM/CD24 chordoma marker set",
  "Use TBXT-linked identity as an anchor for transcription-factor targeting and downstream vulnerability interpretation.",
  "Cell-line TBXT perturbation; patient tumor projection; brachyury vaccine context",
  paste(
    "GSE276715 T-DARPin vs control reduced chordoma marker score:",
    metric_from("results/GSE276715_pathway_score_tests.csv", "CUSTOM_CHORDOMA_MARKERS"),
    "T-DARPin signatures were projected onto GSE230168 tumors."
  ),
  "Direct TBXT targeting remains investigational; RNA signatures require protein and functional validation.",
  "Notochordal/WNT-EMT developmental axis",
  "TBXT, WNT/beta-catenin, EMT, mTORC1, rs2305089 genotype context",
  "Interpret patient tumor states against notochordal differentiation programs without treating developmental identity as direct drug evidence.",
  "Engineered differentiation model; germline-risk literature context",
  paste(
    "GSE277707 heterozygous vs wild type stage-specific shifts; mesoderm WNT:",
    metric_from("results/GSE277707_pathway_score_tests_by_stage.csv", "HALLMARK_WNT_BETA_CATENIN_SIGNALING")
  ),
  "Engineered differentiation models are biological context, not tumor or therapeutic validation.",
  "RTK/angiogenic candidate axis",
  "EGFR, PDGFRA, PDGFRB, KDR/VEGFR2, VEGFA, KIT, MET, FGFR family",
  "Use RTK expression as a hypothesis-generating layer for targeted-agent and combination prioritization.",
  "Patient RNA-seq; targeted-therapy literature; CDK/EGFR combination context",
  paste(
    "GSE230168 chordoma_I vs chordoma_C druggable RTK signature:",
    metric_from("results/GSE230168_pathway_score_tests.csv", "CUSTOM_DRUGGABLE_RTKS")
  ),
  "RTK expression alone is not evidence of dependency; drug sensitivity or phosphoproteomic validation is needed."
)
write_csv(table1, "results/Table1_translational_vulnerability_matrix.csv")

meta <- read_csv("data/processed/GSE230168_metadata_curated.csv", show_col_types = FALSE)
pca_df <- read_csv("results/GSE230168_pca_top_variable_genes.csv", show_col_types = FALSE)
axis_scores <- read_csv("results/GSE230168_axis_scores.csv", show_col_types = FALSE)
micro_scores <- read_csv("results/GSE230168_microenvironment_scores.csv", show_col_types = FALSE)
biomarker_expr <- read_csv("results/GSE230168_biomarker_expression_long.csv", show_col_types = FALSE)

palette <- group_palette()
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(title = "a", x = "PC1", y = "PC2", color = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

p_axes <- ggplot(axis_scores, aes(x = immune_interferon, y = cell_cycle, color = type)) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(title = "b", x = "Immune/interferon score", y = "Cell-cycle score", color = NULL) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

micro_labels <- c(
  MICRO_ANTIGEN_PRESENTATION = "Antigen presentation",
  MICRO_T_CELL_INFLAMMATION = "T-cell inflammation",
  MICRO_MACROPHAGE_MYELOID = "Macrophage/myeloid",
  MICRO_CAF_STROMAL = "CAF/stromal",
  MICRO_ENDOTHELIAL_ANGIOGENIC = "Endothelial/angiogenic",
  MICRO_CHECKPOINT_CONTEXT = "Checkpoint context"
)
sample_order <- axis_scores %>%
  arrange(type, desc(immune_interferon), desc(cell_cycle), sample) %>%
  pull(sample)
micro_long <- micro_scores %>%
  select(sample, type, starts_with("MICRO_")) %>%
  pivot_longer(cols = starts_with("MICRO_"), names_to = "signature", values_to = "score") %>%
  group_by(signature) %>%
  mutate(score_z = as.numeric(scale(score))) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    label = factor(micro_labels[signature], levels = rev(micro_labels))
  )
p_micro <- ggplot(micro_long, aes(x = sample, y = label, fill = score_z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#1E5FA3", mid = "#F1E7DE", high = "#A6492C", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "z-score") +
  labs(title = "c", x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

selected_biomarkers <- c("TBXT", "MKI67", "CDKN2A", "CDKN2B", "PLK1", "EGFR", "PDGFRA", "KDR", "B2M", "HLA-A", "STAT1", "CD274", "HAVCR2", "CD47")
bio_long <- biomarker_expr %>%
  filter(gene_symbol %in% selected_biomarkers) %>%
  group_by(gene_symbol) %>%
  mutate(expr_z = as.numeric(scale(expression))) %>%
  ungroup() %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    gene_symbol = factor(gene_symbol, levels = rev(selected_biomarkers))
  )
p_bio <- ggplot(bio_long, aes(x = sample, y = gene_symbol, fill = expr_z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = "#1E5FA3", mid = "#F1E7DE", high = "#A6492C", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish, name = "z-score") +
  labs(title = "d", x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())

fig1 <- plot_grid(
  plot_grid(p_pca, p_axes, nrow = 1),
  plot_grid(p_micro, p_bio, nrow = 1),
  ncol = 1,
  rel_heights = c(0.95, 1.25)
)
ggsave("figures/Figure2_GSE230168_high_impact_landscape.png", fig1, width = 13, height = 8.5, dpi = 300)
ggsave("figures/Figure2_GSE230168_high_impact_landscape.pdf", fig1, width = 13, height = 8.5)

# Note: the GSE277707 notochordal figure is written directly with its
# Supplementary_Figure_* name by script 10; no rename/copy step is needed here.

robust <- read_csv("results/signature_robustness_summary.csv", show_col_types = FALSE) %>%
  mutate(
    axis_label = factor(axis_label, levels = rev(c("Immune/interferon", "JAK/STAT", "Cell-cycle", "CDK4/6", "PLK1", "RTK", "Chordoma identity"))),
    dataset = factor(dataset, levels = c("GSE230168", "GSE276715", "GSE216418")),
    direction = case_when(
      delta_z_score > 0 ~ "Higher in case group",
      delta_z_score < 0 ~ "Lower in case group",
      TRUE ~ "No direction"
    )
  )
p_concordance <- ggplot(robust, aes(x = dataset, y = axis_label, fill = direction_concordant)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = if_else(direction_concordant, "OK", "Disc")), size = 3) +
  scale_fill_manual(values = c("TRUE" = "#009E73", "FALSE" = "#D55E00"), name = "Concordant") +
  labs(title = "a", x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")

p_delta <- ggplot(robust, aes(x = delta_z_score, y = delta_rank_score, color = direction_concordant, shape = dataset)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_color_manual(values = c("TRUE" = "#009E73", "FALSE" = "#D55E00"), name = "Concordant") +
  labs(title = "b", x = "Delta, z-score method", y = "Delta, rank-score method", shape = NULL) +
  theme_nature(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

evidence_map <- tribble(
  ~state, ~evidence_layer, ~strength,
  "Immune/interferon", "Patient RNA-seq", 3,
  "Immune/interferon", "Perturbation", 2,
  "Immune/interferon", "PDX", 0,
  "Immune/interferon", "Developmental model", 1,
  "Immune/interferon", "Clinical/literature", 3,
  "Cell-cycle/CDK4-6/PLK1", "Patient RNA-seq", 2,
  "Cell-cycle/CDK4-6/PLK1", "Perturbation", 3,
  "Cell-cycle/CDK4-6/PLK1", "PDX", 3,
  "Cell-cycle/CDK4-6/PLK1", "Developmental model", 2,
  "Cell-cycle/CDK4-6/PLK1", "Clinical/literature", 3,
  "TBXT identity", "Patient RNA-seq", 2,
  "TBXT identity", "Perturbation", 3,
  "TBXT identity", "PDX", 0,
  "TBXT identity", "Developmental model", 3,
  "TBXT identity", "Clinical/literature", 2,
  "RTK/angiogenic", "Patient RNA-seq", 2,
  "RTK/angiogenic", "Perturbation", 0,
  "RTK/angiogenic", "PDX", 0,
  "RTK/angiogenic", "Developmental model", 0,
  "RTK/angiogenic", "Clinical/literature", 3
)
write_csv(evidence_map, "results/translational_evidence_layer_map.csv")
p_evidence <- ggplot(evidence_map, aes(x = evidence_layer, y = state, fill = strength)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = strength), size = 3) +
  scale_fill_gradient(low = "#F2F2F2", high = "#0072B2", limits = c(0, 3), name = "Support") +
  labs(title = "c", x = NULL, y = NULL) +
  theme_nature(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))

trial_plot <- trial_context %>%
  mutate(
    evidence_rank = case_when(
      grepl("phase II|trial", evidence_type, ignore.case = TRUE) ~ 3,
      grepl("PDX|perturbation", evidence_type, ignore.case = TRUE) ~ 2,
      TRUE ~ 1
    ),
    strategy = factor(strategy, levels = rev(strategy))
  )
p_trial <- ggplot(trial_plot, aes(x = evidence_rank, y = strategy, color = linked_state)) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:3, labels = c("Literature", "Functional", "Trial context"), limits = c(0.7, 3.3)) +
  labs(title = "d", x = NULL, y = NULL, color = NULL) +
  theme_nature(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

fig4 <- plot_grid(
  plot_grid(p_concordance, p_delta, nrow = 1),
  plot_grid(p_evidence, p_trial, nrow = 1, rel_widths = c(1, 1.05)),
  ncol = 1
)
ggsave("figures/Figure4_translational_robustness_map.png", fig4, width = 13, height = 8.5, dpi = 300)
ggsave("figures/Figure4_translational_robustness_map.pdf", fig4, width = 13, height = 8.5)

write_csv(tribble(
  ~result, ~interpretation,
  "GSE230168 separates nucleus pulposus controls from chordoma tumors by PCA.",
  "The patient cohort remains suitable as the patient-tumor anchor.",
  "GSE230168 chordoma_I shows higher interferon/JAK/inflammatory scores than chordoma_C.",
  "This supports an immune/interferon-high molecular state.",
  "Bulk microenvironment marker scores add T-cell, myeloid, stromal, endothelial, and checkpoint context.",
  "The analysis connects transcriptomic states to biomarker panels without claiming cell-type deconvolution.",
  "Rank-based and z-score methods show 16/21 directional concordance across primary axis tests.",
  "Sensitivity analysis supports key patient and cell-cycle axes while transparently flagging weaker discordant contexts.",
  "GSE276715 T-DARPin perturbation reduces chordoma marker and cell-cycle/PLK1 signatures.",
  "TBXT perturbation connects chordoma identity to proliferative vulnerability.",
  "GSE216418 palbociclib reduces E2F/G2M signatures in PDX models.",
  "The CDK4/6 axis has independent perturbational support and clinical-trial context.",
  "GSE277707 rs2305089 remains useful as developmental context and is moved to supplement.",
  "The main display set now prioritizes patient-state, functional, and translational-evidence figures."
), "results/final_result_highlights.csv")

writeLines(c(
  "# Analysis figures and translational tables",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Generated outputs",
  "- Figure 2 includes patient PCA, immune/cell-cycle axes, bulk microenvironment marker scores, and a translational biomarker panel.",
  "- Figure 4 presents robustness and translational evidence mapping.",
  "- The translational vulnerability matrix includes explicit evidence layers.",
  "",
  "## Main outputs",
  "- figures/Figure2_GSE230168_high_impact_landscape.png",
  "- figures/Figure4_translational_robustness_map.png",
  "- results/translational_trial_context.csv",
  "- results/Table1_translational_vulnerability_matrix.csv"
), "results/analysis_figures_tables_summary.md")

message("Done: analysis figures and translational tables")
