#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Generate validation and quality-control outputs for the analysis.
#
# Main inputs:
#   - data/raw/Yin2024_mmc2_clinical_cluster.xlsx
#   - results/GSE230168_axis_scores.csv
#   - results/GSE239531_axis_scores.csv
#   - results/Yin2024_unsup_axis_scores.csv
#   - previously generated score/test tables in results/
#
# Main outputs:
#   - results/yin_survival_reconciliation.csv
#   - results/axis_correlation_bootstrap_ci.csv
#   - results/tested_signature_inventory.csv
#   - results/custom_signature_genes.csv
#   - results/label_transfer_qc.csv
#   - results/proteomic_label_transfer_recovery_qc.csv
#   - results/proteomic_label_transfer_recovery_addendum.md
#   - results/validation_qc_summary.md
#
# Statistical approach:
#   Reconciles published Yin survival summaries, computes Spearman correlations
#   with bootstrap confidence intervals, inventories tested signatures, and
#   exports custom gene sets, label-transfer QC, and Yin recovery-audit outputs
#   for transparent reporting.
#
# Interpretation notes:
#   This script documents limitations explicitly: weak correlations, small
#   cohorts, mixed rank-score sensitivity, unavailable Yin channel-patient key,
#   and hypothesis-generating claims.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

suppressPackageStartupMessages({
  library(readxl)
})

set.seed(1801)

read_processed_scores <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }
  read_csv(path, show_col_types = FALSE)
}

safe_round <- function(x, digits = 3) {
  ifelse(is.na(x), NA_real_, round(x, digits))
}

# ---- 1. Yin clinical-table survival reconciliation -------------------
yin_clin <- read_excel("data/raw/Yin2024_mmc2_clinical_cluster.xlsx", sheet = 1)
yin_cluster_labels <- c("1" = "BME", "2" = "MES", "3" = "MET")
yin_cluster_names <- c(
  "1" = "Bone-microenvironment-dominant / immunogenic",
  "2" = "Mesenchymal-derived",
  "3" = "MET-mediated / RTK-linked"
)

survival_recon <- yin_clin %>%
  mutate(
    cluster_id = as.character(.data[["Cluster"]]),
    subtype = recode(cluster_id, !!!yin_cluster_labels),
    subtype_description = recode(cluster_id, !!!yin_cluster_names),
    survival_months = as.numeric(.data[["survival time"]]),
    overall_survival_raw = as.character(.data[["Overall survival"]]),
    death_event = !grepl("alive", tolower(overall_survival_raw))
  ) %>%
  group_by(cluster_id, subtype, subtype_description) %>%
  summarise(
    n_patients = n(),
    n_with_survival_time = sum(!is.na(survival_months)),
    mean_survival_months = mean(survival_months, na.rm = TRUE),
    median_survival_months = median(survival_months, na.rm = TRUE),
    min_survival_months = min(survival_months, na.rm = TRUE),
    max_survival_months = max(survival_months, na.rm = TRUE),
    death_events = sum(death_event, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    authoritative_source = "Yin2024_mmc2 clinical table; per-patient published cluster labels",
    reporting_note = if_else(
      subtype == "BME",
      "Use 61.4 months as the mean survival for the immunogenic/BME subtype; remove the previous incorrect mean-survival claim.",
      "Report as comparator subtype from the same published clinical table."
    )
  )
write_csv(survival_recon, "results/yin_survival_reconciliation.csv")

# ---- 2. Bootstrap CIs for weak immune-cell-cycle correlations --------
bootstrap_spearman <- function(df, dataset, cohort_note, x_col = "immune_interferon",
                               y_col = "cell_cycle", n_boot = 2000) {
  sub <- df %>%
    select(all_of(c(x_col, y_col))) %>%
    filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))

  observed <- suppressWarnings(cor(sub[[x_col]], sub[[y_col]], method = "spearman"))
  p_value <- suppressWarnings(cor.test(sub[[x_col]], sub[[y_col]],
                                       method = "spearman", exact = FALSE)$p.value)
  boot <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(sub)), replace = TRUE)
    suppressWarnings(cor(sub[[x_col]][idx], sub[[y_col]][idx], method = "spearman"))
  })
  boot <- boot[is.finite(boot)]

  tibble(
    dataset = dataset,
    cohort_note = cohort_note,
    x_axis = x_col,
    y_axis = y_col,
    n = nrow(sub),
    spearman_rho = observed,
    p_value = p_value,
    bootstrap_resamples = n_boot,
    ci_lower = unname(quantile(boot, 0.025, na.rm = TRUE)),
    ci_upper = unname(quantile(boot, 0.975, na.rm = TRUE)),
    interpretation = case_when(
      abs(observed) < 0.30 ~ "weak directionally consistent correlation; do not call quantitative replication",
      abs(observed) < 0.50 ~ "moderate correlation",
      TRUE ~ "strong correlation"
    )
  )
}

gse230168_axis <- read_processed_scores("results/GSE230168_axis_scores.csv") %>%
  filter(type %in% c("chordoma_I", "chordoma_C"))
gse239531_axis <- read_processed_scores("results/GSE239531_axis_scores.csv")
yin_axis <- read_processed_scores("results/Yin2024_unsup_axis_scores.csv")

axis_cor_ci <- bind_rows(
  bootstrap_spearman(gse230168_axis, "GSE230168", "Patient RNA-seq tumors only"),
  bootstrap_spearman(gse239531_axis, "GSE239531", "Independent patient RNA-seq; immune split is 18/2 and descriptive"),
  bootstrap_spearman(yin_axis, "Yin2024_protein", "Orthogonal clinical-proteomic cohort; TMT channels without patient key")
) %>%
  mutate(across(c(spearman_rho, p_value, ci_lower, ci_upper), ~safe_round(.x, 3)))
write_csv(axis_cor_ci, "results/axis_correlation_bootstrap_ci.csv")

# ---- 3. Inventory of tested signatures and contrasts -----------------
count_csv_rows <- function(path, dataset, family, contrast_scope) {
  if (!file.exists(path)) {
    return(tibble(
      dataset = dataset, family = family, contrast_scope = contrast_scope,
      source_file = path, n_tests = NA_integer_, status = "missing"
    ))
  }
  x <- read_csv(path, show_col_types = FALSE)
  tibble(
    dataset = dataset,
    family = family,
    contrast_scope = contrast_scope,
    source_file = path,
    n_tests = nrow(x),
    status = "available"
  )
}

inventory <- bind_rows(
  count_csv_rows("results/GSE230168_pathway_score_tests.csv", "GSE230168", "selected Hallmark/custom pathways", "chordoma_I vs chordoma_C and tumor vs nucleus pulposus"),
  count_csv_rows("results/GSE230168_microenvironment_score_tests.csv", "GSE230168", "bulk microenvironment marker scores", "chordoma_I vs chordoma_C"),
  count_csv_rows("results/GSE230168_biomarker_group_tests.csv", "GSE230168", "single-gene translational biomarker panel", "chordoma_I vs chordoma_C"),
  count_csv_rows("results/GSE276715_pathway_score_tests.csv", "GSE276715", "selected Hallmark/custom pathways", "T-DARPin vs control"),
  count_csv_rows("results/GSE216418_pathway_score_tests.csv", "GSE216418", "selected Hallmark/custom pathways", "palbociclib vs untreated"),
  count_csv_rows("results/GSE277707_pathway_score_tests_all.csv", "GSE277707", "selected Hallmark/custom pathways", "rs2305089 heterozygous vs wild type"),
  count_csv_rows("results/GSE239531_axis_group_tests.csv", "GSE239531", "axis scores", "data-driven immune_high vs immune_low; descriptive 18/2 split"),
  count_csv_rows("results/GSE239531_pathway_score_tests.csv", "GSE239531", "selected Hallmark/custom pathways", "data-driven immune_high vs immune_low; descriptive 18/2 split"),
  count_csv_rows("results/Yin2024_supervised_axis_kruskal.csv", "Yin2024_protein", "axis scores", "marker-label transferred subtype Kruskal-Wallis"),
  count_csv_rows("results/Yin2024_supervised_pathway_kruskal.csv", "Yin2024_protein", "selected Hallmark/custom pathways", "marker-label transferred subtype Kruskal-Wallis"),
  count_csv_rows("results/signature_robustness_summary.csv", "multi-dataset", "z-score vs rank-based score sensitivity", "21 primary axis-direction tests"),
  count_csv_rows("results/effect_size_stability_summary.csv", "multi-dataset", "Cliff's delta/bootstrap/leave-one-out stability", "20 prespecified core contrasts")
) %>%
  mutate(
    reporting_role = case_when(
      grepl("18/2", contrast_scope) ~ "descriptive context only",
      grepl("21 primary", contrast_scope) ~ "secondary sensitivity output; global robustness threshold not met",
      TRUE ~ "primary or secondary analysis as labeled"
    )
  )
write_csv(inventory, "results/tested_signature_inventory.csv")

# ---- 4. Full custom signature gene table -----------------------------
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

custom_sets <- selected_gene_sets()[grepl("^CUSTOM_", names(selected_gene_sets()))]
micro_sets <- microenvironment_sets()
sym_map <- read_csv("data/processed/symbol2entrez.csv", show_col_types = FALSE) %>%
  transmute(gene_symbol = as.character(symbol), entrez_id = as.character(entrez_id)) %>%
  filter(!is.na(entrez_id), entrez_id != "") %>%
  distinct(entrez_id, .keep_all = TRUE)

table_from_sets <- function(sets, family, axis_or_context, source_note) {
  bind_rows(lapply(names(sets), function(set_name) {
    tibble(
      family = family,
      signature_name = set_name,
      axis_or_context = axis_or_context[set_name],
      entrez_id = as.character(sets[[set_name]]),
      source_note = source_note
    )
  }))
}

custom_axis_context <- c(
  CUSTOM_ANTIGEN_PRESENTATION = "antigen presentation / immune state",
  CUSTOM_CHORDOMA_MARKERS = "chordoma identity",
  CUSTOM_CDK46_AXIS = "cell-cycle / CDK4-6",
  CUSTOM_PLK1_AXIS = "cell-cycle / PLK1",
  CUSTOM_DRUGGABLE_RTKS = "RTK/angiogenic context"
)
micro_context <- c(
  MICRO_ANTIGEN_PRESENTATION = "bulk antigen-presentation marker context",
  MICRO_T_CELL_INFLAMMATION = "bulk T-cell inflammation marker context",
  MICRO_MACROPHAGE_MYELOID = "bulk macrophage/myeloid marker context",
  MICRO_CAF_STROMAL = "bulk CAF/stromal marker context",
  MICRO_ENDOTHELIAL_ANGIOGENIC = "bulk endothelial/angiogenic marker context",
  MICRO_CHECKPOINT_CONTEXT = "bulk checkpoint marker context"
)

custom_gene_table <- bind_rows(
  table_from_sets(custom_sets, "custom pathway/signature", custom_axis_context, "Defined in scripts/00_chordoma_utils.R"),
  table_from_sets(micro_sets, "bulk microenvironment marker score", micro_context, "Defined in scripts/02_gse230168_microenvironment_biomarker_context.R"),
  biomarker_panel %>%
    transmute(
      family = "single-gene translational biomarker panel",
      signature_name = paste0("BIOMARKER_", gene_symbol),
      axis_or_context = axis,
      entrez_id = as.character(entrez_id),
      source_note = "Defined in scripts/02_gse230168_microenvironment_biomarker_context.R"
    )
) %>%
  left_join(sym_map, by = "entrez_id") %>%
  mutate(
    gene_symbol = coalesce(gene_symbol, paste0("Entrez:", entrez_id)),
    scoring_note = case_when(
      family == "single-gene translational biomarker panel" ~ "single-gene expression, not gene-set score",
      grepl("^MICRO_", signature_name) ~ "bulk marker-set score; not cell-type deconvolution",
      TRUE ~ "mean sample-wise z-score; rank-based sensitivity where applicable"
    )
  ) %>%
  select(family, signature_name, axis_or_context, gene_symbol, entrez_id, source_note, scoring_note) %>%
  arrange(family, signature_name, gene_symbol)
write_csv(custom_gene_table, "results/custom_signature_genes.csv")

# ---- 5. Yin label-transfer QC ----------------------------------------
assigned_counts <- read_processed_scores("results/Yin2024_channel_cluster_assignment.csv") %>%
  count(cluster, name = "assigned_channels") %>%
  mutate(cluster_id = recode(cluster, BME = "1", MES = "2", MET = "3"))

published_counts <- yin_clin %>%
  mutate(cluster_id = as.character(.data[["Cluster"]])) %>%
  count(cluster_id, name = "published_patients") %>%
  mutate(cluster = recode(cluster_id, !!!yin_cluster_labels))

label_transfer_qc <- full_join(assigned_counts, published_counts, by = c("cluster_id", "cluster")) %>%
  mutate(
    subtype_description = recode(cluster_id, !!!yin_cluster_names),
    count_difference = assigned_channels - published_patients,
    percent_difference_vs_published = 100 * count_difference / published_patients,
    assignment_rule = "nearest subtype-marker centroid using Yin mmc4 markers z-scored across TMT channels",
    marker_overlap_with_tested_signatures = "none for immune/cell-cycle/RTK/chordoma identity test signatures by construction",
    verifiable_per_sample_accuracy_available = FALSE,
    limitation = "No channel-to-patient key is published, so label-transfer accuracy cannot be verified per patient; survival uses published clinical labels only."
  ) %>%
  mutate(across(c(percent_difference_vs_published), ~safe_round(.x, 1))) %>%
  select(
    cluster_id, cluster, subtype_description, assigned_channels, published_patients,
    count_difference, percent_difference_vs_published, assignment_rule,
    marker_overlap_with_tested_signatures, verifiable_per_sample_accuracy_available,
    limitation
  )
write_csv(label_transfer_qc, "results/label_transfer_qc.csv")

# ---- 6. Yin label-transfer recovery audit ----------------------------
protein_header <- read_csv("data/raw/Yin2024_protein_abundance.csv", n_max = 0, show_col_types = FALSE)
xml_files <- c(
  "data/raw/IPX0009051000/IPX0009051000.xml",
  "data/raw/IPX0009051000/PX_IPX0009051000.xml"
)
xml_available <- xml_files[file.exists(xml_files)]
protein_channel_examples <- paste(head(setdiff(names(protein_header), names(protein_header)[1]), 8), collapse = "; ")
clinical_sample_examples <- paste(head(as.character(yin_clin[["Sample"]]), 5), collapse = "; ")

recovery_qc <- tibble(
  source_checked = c(
    "iProX XML metadata",
    "protein abundance header",
    "Yin mmc2 clinical table",
    "Yin mmc4 marker centroids",
    "local raw archive files"
  ),
  local_path_or_pattern = c(
    ifelse(length(xml_available) > 0, paste(xml_available, collapse = "; "), paste(xml_files, collapse = "; ")),
    "data/raw/Yin2024_protein_abundance.csv",
    "data/raw/Yin2024_mmc2_clinical_cluster.xlsx",
    "data/raw/Yin2024_mmc4_cluster_markers.xlsx",
    "data/raw/IPX0009051000/IPX0009051001/*.partial; data/raw/IPX0009051000/IPX0009051001/*.aspera-ckpt"
  ),
  key_signal = c(
    "dataset files and iProX URLs only; no TMT channel-to-patient key",
    paste0("TMT channel headers, for example ", protein_channel_examples),
    paste0("patient samples, for example ", clinical_sample_examples, "; no TMT channel columns"),
    "45 subtype marker centroids only",
    "partial/checkpoint files only; no completed archive metadata available locally"
  ),
  recovered_key = FALSE
)
write_csv(recovery_qc, "results/proteomic_label_transfer_recovery_qc.csv")

writeLines(c(
  "# Proteomic label-transfer recovery addendum",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Sources re-inspected",
  "- iProX XML metadata under data/raw/IPX0009051000/.",
  "- Yin protein abundance header.",
  "- Yin mmc2 clinical table.",
  "- Yin mmc4 subtype-marker centroids.",
  "- Local iProX partial/checkpoint files.",
  "",
  "## Result",
  "- recovered_key = FALSE.",
  "- The local original files do not contain a channel-to-patient key linking TMT channels (for example F1_127N) to clinical samples (for example SA1).",
  "- Supervised proteomic subtype tests therefore remain marker-label-transfer analyses and should be interpreted directionally.",
  "- Yin survival statements use the published per-patient clinical table, not transferred TMT-channel labels."
), "results/proteomic_label_transfer_recovery_addendum.md")

# ---- 7. Summary ------------------------------------------------------
writeLines(c(
  "# Validation and quality-control outputs",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Validation checks",
  "- Yin survival is reconciled to the published clinical table: BME mean 61.4 months, MES 44.9 months, MET 42.6 months.",
  "- Immune-cell-cycle correlations are weak and directionally consistent; bootstrap CIs are in results/axis_correlation_bootstrap_ci.csv.",
  "- The 16/21 rank-score concordance result remains below the prespecified 80% threshold and should be treated as mixed sensitivity evidence.",
  "- Custom gene sets, microenvironment marker sets, and biomarker panel genes are exported in full.",
  "- Yin label transfer is approximate; no per-sample accuracy can be verified without a channel-to-patient key.",
  "- A local Yin recovery audit re-inspected iProX XML metadata, protein headers, clinical clusters, and marker centroids; recovered_key = FALSE.",
  "",
  "## Outputs",
  "- results/major_revision_survival_reconciliation.csv",
  "- results/axis_correlation_bootstrap_ci.csv",
  "- results/tested_signature_inventory.csv",
  "- results/custom_signature_genes.csv",
  "- results/label_transfer_qc.csv",
  "- results/proteomic_label_transfer_recovery_qc.csv",
  "- results/proteomic_label_transfer_recovery_addendum.md"
), "results/validation_qc_summary.md")

message("Done: validation and quality-control outputs")
