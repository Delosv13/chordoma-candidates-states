#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Audit Yin 2024 proteomic inputs before confirmatory context checks.
#
# Main inputs:
#   - data/raw/Yin2024_protein_abundance.csv
#   - data/raw/Yin2024_mmc2_clinical_cluster.xlsx
#   - data/raw/Yin2024_mmc4_cluster_markers.xlsx
#
# Main outputs:
#   - results/Yin2024_published_cluster_sizes.csv
#   - results/Yin2024_input_qc.csv
#   - results/Yin2024_input_qc_summary.md
#
# Statistical approach:
#   Input audit only. Counts TMT sample channels, published clinical clusters,
#   marker-centroid rows, and records that no channel-to-patient key is public.
#
# Interpretation notes:
#   This script explains why Yin context is split: script 08 is unsupervised,
#   and script 09 uses non-overlapping marker centroids for label transfer.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

suppressPackageStartupMessages({
  library(readxl)
})

required_files <- c(
  protein_abundance = "data/raw/Yin2024_protein_abundance.csv",
  clinical_table = "data/raw/Yin2024_mmc2_clinical_cluster.xlsx",
  marker_table = "data/raw/Yin2024_mmc4_cluster_markers.xlsx"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing Yin 2024 input(s): ", paste(missing_files, collapse = ", "))
}

prot_header <- names(read_csv(required_files[["protein_abundance"]], n_max = 0, show_col_types = FALSE))
sample_cols <- setdiff(prot_header, prot_header[1])
sample_cols <- sample_cols[!grepl("RefInt", sample_cols, ignore.case = TRUE)]

clinical <- read_excel(required_files[["clinical_table"]], sheet = 1)
cluster_summary <- clinical %>%
  mutate(
    cluster_id = as.character(.data[["Cluster"]]),
    subtype = recode(cluster_id, "1" = "BME", "2" = "MES", "3" = "MET")
  ) %>%
  count(cluster_id, subtype, name = "published_patients") %>%
  arrange(cluster_id)

marker_tbl <- read_excel(required_files[["marker_table"]], sheet = 1)

input_qc <- tibble(
  item = c(
    "protein_sample_tmt_channels",
    "clinical_patients",
    "subtype_marker_rows",
    "channel_to_patient_key_available"
  ),
  value = c(
    as.character(length(sample_cols)),
    as.character(nrow(clinical)),
    as.character(nrow(marker_tbl)),
    "FALSE"
  ),
  note = c(
    "Reference-pool columns excluded by script 08.",
    "Published per-patient clinical table used for survival reconciliation.",
    "Yin subtype-marker centroids used by script 09 for non-circular label transfer.",
    "No channel-to-patient key is published; per-sample label-transfer accuracy cannot be verified."
  )
)

write_csv(cluster_summary, "results/Yin2024_published_cluster_sizes.csv")
write_csv(input_qc, "results/Yin2024_input_qc.csv")

writeLines(c(
  "# Yin 2024 proteomic input QC",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Status",
  "- Required Yin protein abundance, clinical, and marker files are present.",
  paste0("- Sample TMT channels available for unsupervised scoring: ", length(sample_cols), "."),
  paste0("- Published clinical patients: ", nrow(clinical), "."),
  "- No channel-to-patient key is published, so survival must use the published clinical table rather than label-transferred TMT channels.",
  "",
  "## Next scripts",
  "- Script 15 scores Yin protein data without subtype labels.",
  "- Script 16 performs marker-based subtype label transfer for directional cross-modality testing.",
  "- Script 18 reconciles published survival values and label-transfer QC for the major revision."
), "results/Yin2024_input_qc_summary.md")

message("Done: Yin 2024 input QC")
