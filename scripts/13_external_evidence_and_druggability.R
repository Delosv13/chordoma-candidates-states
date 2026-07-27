#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Align transcriptomic states with external literature and druggability.
#
# Main inputs:
#   - results/bootstrap_pathway_delta_ci.csv
#   - results/GSE230168_* contextual summaries when available
#   - curated public literature/trial identifiers encoded in this script
#
# Main outputs:
#   - results/external_evidence_alignment.csv
#   - results/translational_druggability_expanded.csv
#   - results/translational_trial_context.csv
#
# Statistical approach:
#   No new inferential model is fit. The script reports prespecified
#   effect summaries alongside external evidence layers and limitations.
#
# Interpretation notes:
#   External studies are used as alignment context only. They are not reanalyzed
#   and should not be interpreted as independent statistical validation here.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

metric_from <- function(path, feature_id, delta_col = "delta_case_minus_control") {
  if (!file.exists(path)) return("not available")
  x <- read_csv(path, show_col_types = FALSE)
  row <- x %>% filter(.data$feature == feature_id) %>% slice_head(n = 1)
  if (nrow(row) == 0) return("not available")
  paste0(
    "delta=", signif(row[[delta_col]], 3),
    "; 95% CI ",
    signif(row$delta_ci_low, 3),
    " to ",
    signif(row$delta_ci_high, 3),
    "; Cliff's delta=",
    signif(row$cliffs_delta, 3)
  )
}

message("Building external evidence and druggability alignment tables")

external_evidence <- tribble(
  ~state, ~our_axis_or_marker, ~our_result, ~external_evidence_type, ~external_source, ~source_identifier, ~alignment, ~limitation,
  "Immune/interferon-high",
  "IL6-JAK-STAT3; interferon; T-cell inflammation; checkpoint context",
  paste(
    "GSE230168 chordoma_I immune/interferon stability:",
    metric_from("results/bootstrap_pathway_delta_ci.csv", "immune_interferon")
  ),
  "Spatial and immune profiling",
  "Multimodal profiling of chordoma immunity reveals distinct immune contextures",
  "PubMed 41586579",
  "Supports the biological plausibility of immune-context heterogeneity in chordoma.",
  "External immune/spatial studies are not reanalyzed here and are used only as alignment evidence.",
  "Immune/interferon-high",
  "T-cell inflammation; myeloid/macrophage; checkpoint markers",
  "GSE230168 chordoma_I showed higher T-cell inflammation, macrophage/myeloid, and checkpoint-context marker scores.",
  "Methylation, spatial transcriptomics, and multiplexed immunofluorescence",
  "Molecular profiling of sacral chordomas through methylation, spatial transcriptomics and multiplexed immunofluorescence",
  "PubMed 41924303",
  "Supports the need to connect methylation/transcriptional classes with spatial immune and stromal context.",
  "Different anatomic site and external study design limit direct comparability.",
  "Cell-cycle/CDK4-6/PLK1-high",
  "E2F; G2M; CDK4/6; PLK1",
  paste(
    "GSE216418 palbociclib cell-cycle stability:",
    metric_from("results/bootstrap_pathway_delta_ci.csv", "cell_cycle")
  ),
  "PDX perturbation and clinical trial context",
  "Palbociclib in advanced chordoma",
  "GSE216418; NCT03110744",
  "Supports CDK4/6 pathway scoring as a translational hypothesis linked to trial development.",
  "Transcriptomic response signatures are not patient response predictors.",
  "Cell-cycle/CDK4-6/PLK1-high",
  "CDK/EGFR overlap",
  "GSE230168 RTK scores and GSE216418 cell-cycle effects support a combination-prioritization hypothesis.",
  "Preclinical therapeutic strategy",
  "Preclinical CDK/EGFR therapeutic strategy in chordoma",
  "PubMed 41340469",
  "Supports evaluating proliferative and RTK states together rather than as isolated axes.",
  "Preclinical literature alignment does not establish a clinical combination biomarker.",
  "TBXT/chordoma identity",
  "TBXT; chordoma marker score; T-DARPin perturbation",
  paste(
    "GSE276715 chordoma marker stability:",
    metric_from("results/bootstrap_pathway_delta_ci.csv", "CUSTOM_CHORDOMA_MARKERS")
  ),
  "Functional TBXT perturbation and vaccine context",
  "Selective targeting of TBXT with DARPins; yeast-brachyury vaccine study",
  "PubMed 40901948; NCT02383498",
  "Supports TBXT as a lineage and therapeutic-development anchor.",
  "Brachyury-directed strategies remain investigational and require protein/functional validation.",
  "Proteomic/precision-treatment context",
  "Immune, RTK, and cell-cycle state interpretation",
  "The current study uses transcriptomic states and evidence layers.",
  "Clinical proteomic classification",
  "Clinical-proteomic classification and precision treatment strategy of chordoma",
  "PubMed 39368483",
  "Supports multi-layer therapeutic-state framing beyond RNA alone.",
  "Proteomic data are not reanalyzed, so claims are limited to concordant rationale.",
  "Developmental/notochordal context",
  "TBXT rs2305089; WNT/EMT; stage-specific cell-cycle shifts",
  "GSE277707 is retained as developmental context rather than therapeutic validation.",
  "Genetic/developmental risk context",
  "TBXT rs2305089 and chordoma/notochordal tumor literature",
  "PubMed 42121852; GSE277707",
  "Supports separating lineage/developmental identity from druggable tumor stress states.",
  "Engineered differentiation models are not tumors."
)

druggability <- tribble(
  ~state, ~biomarker_or_axis, ~target, ~agent_or_class, ~evidence_source, ~evidence_level, ~our_result, ~recommended_validation, ~limitation,
  "Cell-cycle/CDK4-6/PLK1-high",
  "E2F/G2M/CDK4-6 score; CDKN2A/B; RB pathway context",
  "CDK4/6",
  "Palbociclib or CDK4/6 inhibitor class",
  "GSE216418; NCT03110744",
  "PDX perturbation plus clinical-trial context",
  metric_from("results/bootstrap_pathway_delta_ci.csv", "cell_cycle"),
  "Test whether baseline cell-cycle and RB/CDKN2A-B status associate with response in biomarker-rich cohorts.",
  "No validated transcriptomic response predictor.",
  "Cell-cycle/CDK4-6/PLK1-high",
  "PLK1 axis",
  "PLK1",
  "Volasertib or PLK1 inhibitor class",
  "GSE276715; chordoma PDX literature",
  "Perturbational and preclinical context",
  metric_from("results/bootstrap_pathway_delta_ci.csv", "CUSTOM_PLK1_AXIS"),
  "Validate PLK1 protein/activity and drug sensitivity in models with high PLK1-axis scores.",
  "RNA score is not dependency evidence.",
  "TBXT/chordoma identity",
  "TBXT and chordoma marker score",
  "TBXT/brachyury",
  "T-DARPin strategy; brachyury vaccine context",
  "GSE276715; PubMed 40901948; NCT02383498",
  "Functional perturbation plus immunotherapy-development context",
  metric_from("results/bootstrap_pathway_delta_ci.csv", "CUSTOM_CHORDOMA_MARKERS"),
  "Confirm TBXT protein suppression and downstream cell-cycle effects in additional chordoma models.",
  "Direct TBXT targeting is investigational.",
  "Immune/interferon-high",
  "IFN/JAK/STAT; T-cell inflammation; checkpoint context",
  "JAK/STAT and immune-checkpoint context",
  "JAK-pathway or immunotherapy-combination hypotheses",
  "GSE230168; immune/spatial profiling literature",
  "Patient RNA-seq plus external immune context",
  metric_from("results/bootstrap_pathway_delta_ci.csv", "immune_interferon"),
  "Use spatial/protein assays to distinguish tumor-cell-intrinsic from microenvironmental immune signal.",
  "Bulk RNA-seq marker scores are not deconvolution.",
  "RTK/angiogenic candidate axis",
  "EGFR; PDGFRA; KDR; VEGFA; druggable RTK score",
  "EGFR/PDGFRA/VEGFR axis",
  "EGFR, PDGFR, VEGFR, or combination targeted strategies",
  "GSE230168; imatinib/lapatinib literature; PubMed 41340469",
  "Patient RNA-seq plus targeted-therapy literature",
  metric_from("results/bootstrap_pathway_delta_ci.csv", "CUSTOM_DRUGGABLE_RTKS"),
  "Validate phosphoproteomic activation and drug sensitivity before therapeutic selection.",
  "RTK expression alone does not prove dependency.",
  "Developmental/notochordal context",
  "TBXT rs2305089; WNT/EMT/mTORC1",
  "Developmental lineage programs",
  "Contextual biomarker interpretation",
  "GSE277707; rs2305089 literature",
  "Developmental model context",
  "Stage-specific pathway shifts in GSE277707.",
  "Use as biological context when interpreting tumor states.",
  "Not a therapeutic validation model."
)

write_csv(external_evidence, "results/external_evidence_alignment.csv")
write_csv(druggability, "results/translational_druggability_expanded.csv")
write_csv(druggability, "results/translational_trial_context.csv")

writeLines(c(
  "# External Evidence And Druggability Alignment",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Method",
  "- External spatial, immune, proteomic, clinical-trial, and preclinical evidence was aligned to the candidate molecular states.",
  "- External studies are not reanalyzed; they are used as literature/context alignment only.",
  "- Druggability rows map state -> biomarker/axis -> target -> agent/class -> evidence level -> limitation.",
  "",
  "## Main outputs",
  "- results/external_evidence_alignment.csv",
  "- results/translational_druggability_expanded.csv",
  "- results/external_evidence_alignment.csv",
  "- results/translational_druggability_expanded.csv"
), "results/external_evidence_alignment_summary.md")

message("Done: external evidence and druggability alignment")
