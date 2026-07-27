#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Reproject chordoma state signatures into confirmatory RNA context GSE239531.
#
# Main inputs:
#   - data/raw/GSE239531_raw_counts.tsv.gz (downloaded from GEO if absent)
#   - data/raw/NCBI_gene2ensembl.gz or data/raw/hgnc_complete_set.tsv
#   - data/raw/h.all.v2025.1.Hs.entrez.gmt
#
# Main outputs:
#   - data/processed/GSE239531_log2cpm_entrez.csv
#   - results/GSE239531_pathway_scores.csv
#   - results/GSE239531_axis_scores.csv
#   - results/GSE239531_axis_group_tests.csv
#   - results/GSE239531_pathway_score_tests.csv
#   - figures/Supplementary_Figure_GSE239531_independent_validation.*
#
# Statistical approach:
#   Counts are harmonized to Entrez, filtered, transformed as log2(CPM+1), and
#   scored with the same signatures used in GSE230168. An unsupervised k-means
#   split on the immune/interferon axis defines immune_high versus immune_low;
#   group tests use Wilcoxon/BH-FDR and axis association uses Spearman rho.
#
# Interpretation notes:
#   The 18/2 immune split is descriptive and should be presented as external
#   directional confirmatory context rather than a definitive replication cohort.
# -------------------------------------------------------------------------
source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

suppressPackageStartupMessages(library(readr))

gse <- "GSE239531"
raw_gz <- file.path("data/raw", paste0(gse, "_raw_counts.tsv.gz"))
url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=", gse,
  "&format=file&file=", gse, "%5Fraw%5Fcounts%2Etsv%2Egz"
)

# ---- 1. Download processed counts if not already present -------------
if (!file.exists(raw_gz)) {
  message("Downloading ", gse, " raw counts from GEO ...")
  utils::download.file(url, raw_gz, mode = "wb", quiet = FALSE)
}

message("Reading ", gse, " count matrix")
counts_raw <- read_tsv(raw_gz, show_col_types = FALSE)
gene_col <- colnames(counts_raw)[1]
sample_cols <- colnames(counts_raw)[-1]
message("  ", length(sample_cols), " samples, ", nrow(counts_raw), " gene rows")

# ---- 2. Map gene IDs to Entrez (auto-detect ID type) -----------------
ids <- as.character(counts_raw[[gene_col]])
id_type <- if (any(grepl("^ENSG", ids))) {
  "ensembl"
} else if (all(grepl("^[0-9]+$", ids[!is.na(ids)][seq_len(min(50, sum(!is.na(ids))))]))) {
  "entrez"
} else {
  "symbol"
}
message("  Detected gene ID type: ", id_type)

if (id_type == "ensembl") {
  gene_map <- read_human_gene2ensembl()
  mapped <- map_ensembl_matrix_to_entrez(counts_raw, gene_col, sample_cols, gene_map)
  count_matrix <- matrix_from_gene_table(mapped, "entrez_id", sample_cols)
} else if (id_type == "entrez") {
  mapped <- counts_raw %>%
    rename(entrez_id = all_of(gene_col)) %>%
    mutate(entrez_id = as.character(entrez_id)) %>%
    group_by(entrez_id) %>%
    summarise(across(all_of(sample_cols), ~sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
  count_matrix <- matrix_from_gene_table(mapped, "entrez_id", sample_cols)
} else {
  # symbol -> entrez via HGNC complete set already in data/raw
  hgnc <- read_tsv("data/raw/hgnc_complete_set.tsv", show_col_types = FALSE) %>%
    select(symbol, entrez_id) %>%
    filter(!is.na(entrez_id), entrez_id != "") %>%
    mutate(entrez_id = as.character(entrez_id)) %>%
    distinct(symbol, .keep_all = TRUE)
  mapped <- counts_raw %>%
    rename(symbol = all_of(gene_col)) %>%
    inner_join(hgnc, by = "symbol") %>%
    select(entrez_id, all_of(sample_cols)) %>%
    group_by(entrez_id) %>%
    summarise(across(everything(), ~sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
  count_matrix <- matrix_from_gene_table(mapped, "entrez_id", sample_cols)
}

log_cpm <- log2_cpm(count_matrix, min_cpm = 1, min_samples = 3)
message("Retained ", nrow(log_cpm), " expressed Entrez genes")
write_csv(as.data.frame(log_cpm) %>% rownames_to_column("entrez_id"),
          file.path("data/processed", paste0(gse, "_log2cpm_entrez.csv")))

# ---- 3. Score the SAME gene sets / axes used in GSE230168 ------------
gene_sets <- selected_gene_sets()
scores <- score_gene_sets(log_cpm, gene_sets)
axis_scores <- axis_scores_from_pathways(scores)

write_csv(scores, file.path("results", paste0(gse, "_pathway_scores.csv")))
write_csv(gene_set_overlap(log_cpm, gene_sets),
          file.path("results", paste0(gse, "_pathway_gene_overlap.csv")))

# ---- 4. Unsupervised immune-high vs immune-low split -----------------
# van Oost defined two groups by T-cell infiltration; here we recover an
# immune-high vs immune-low split from the immune/interferon axis (no
# external labels needed) and ask whether the two-axis structure holds.
set.seed(1)
imm <- axis_scores$immune_interferon
km <- kmeans(matrix(imm, ncol = 1), centers = 2, nstart = 25)
hi_cluster <- which.max(tapply(imm, km$cluster, mean))
axis_scores$immune_group <- ifelse(km$cluster == hi_cluster, "immune_high", "immune_low")
write_csv(axis_scores, file.path("results", paste0(gse, "_axis_scores.csv")))

meta_split <- axis_scores %>% select(sample, immune_group)

# Replication test 1: does cell-cycle differ between immune groups, and
# what is the immune<->cell-cycle correlation? (compare with GSE230168)
axis_stats <- compare_scores(
  axis_scores %>% select(sample, immune_interferon, cell_cycle, chordoma_marker),
  meta_split, "immune_group", "immune_high", "immune_low"
)
cor_immune_cc <- suppressWarnings(
  cor(axis_scores$immune_interferon, axis_scores$cell_cycle,
      method = "spearman", use = "complete.obs")
)
write_csv(axis_stats, file.path("results", paste0(gse, "_axis_group_tests.csv")))

# Replication test 2: pathway-level contrast immune_high vs immune_low
pathway_stats <- compare_scores(scores, meta_split, "immune_group",
                                "immune_high", "immune_low")
write_csv(pathway_stats, file.path("results", paste0(gse, "_pathway_score_tests.csv")))

# ---- 5. Sanity check: canonical immune markers track the axis --------
immune_markers <- tibble(
  gene_symbol = c("CD274", "PDCD1", "HAVCR2", "STAT1", "CXCL9", "GZMB", "CD8A"),
  entrez_id   = c("29126", "5133", "84868", "6772", "4283", "3002", "925")
)
marker_long <- lapply(seq_len(nrow(immune_markers)), function(i) {
  gid <- immune_markers$entrez_id[i]
  vals <- if (gid %in% rownames(log_cpm)) log_cpm[gid, axis_scores$sample] else rep(NA_real_, nrow(axis_scores))
  tibble(sample = axis_scores$sample, gene_symbol = immune_markers$gene_symbol[i],
         entrez_id = gid, log2cpm = as.numeric(vals))
}) %>% bind_rows() %>%
  left_join(meta_split, by = "sample")
write_csv(marker_long, file.path("results", paste0(gse, "_immune_marker_expression.csv")))

# ---- 6. Figure: external-cohort two-axis landscape -------------------
p_axes <- ggplot(axis_scores, aes(x = immune_interferon, y = cell_cycle, color = immune_group)) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c(immune_high = "#009E73", immune_low = "#D55E00")) +
  labs(title = NULL,
       subtitle = NULL,
       x = "Immune/interferon score", y = "Cell-cycle score", color = NULL) +
  theme_nature(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path("figures", paste0("Supplementary_Figure_", gse, "_independent_validation.png")),
       p_axes, width = 6.5, height = 5, dpi = 300)
ggsave(file.path("figures", paste0("Supplementary_Figure_", gse, "_independent_validation.pdf")),
       p_axes, width = 6.5, height = 5)

# ---- 7. Summary ------------------------------------------------------
n_hi <- sum(axis_scores$immune_group == "immune_high")
imm_axis_row <- axis_stats %>% filter(pathway == "immune_interferon")
writeLines(c(
  paste0("# ", gse, " independent validation (Day 1)"),
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## What this tests",
  "Reprojection of the GSE230168 immune/interferon and cell-cycle axes onto an",
  "independent chordoma RNA-seq cohort (van Oost et al., JITC 2024).",
  "",
  "## Key numbers",
  paste0("- Samples: ", nrow(axis_scores), " (immune_high = ", n_hi,
         ", immune_low = ", nrow(axis_scores) - n_hi, ")."),
  paste0("- Expressed Entrez genes: ", nrow(log_cpm), "."),
  paste0("- Spearman correlation immune vs cell-cycle: ", round(cor_immune_cc, 3), "."),
  if (nrow(imm_axis_row) == 1) paste0(
    "- Immune axis high vs low: delta = ", round(imm_axis_row$delta_a_minus_b, 3),
    ", FDR = ", signif(imm_axis_row$fdr, 3), "."),
  "",
  "## Interpretation guide",
  "- A clear immune_high subgroup + canonical markers (CD274/STAT1/CXCL9/GZMB)",
  "  tracking the immune axis = the immune/interferon state replicates.",
  "- Compare the immune<->cell-cycle correlation sign with GSE230168.",
  "- NOTE: groups here are data-driven, not van Oost's published T-cell labels.",
  "  To run the supervised test, add their infiltration labels (paper Suppl.)",
  "  as a column and rerun compare_scores() on that column.",
  "",
  "## Main outputs",
  paste0("- results/", gse, "_axis_scores.csv"),
  paste0("- results/", gse, "_axis_group_tests.csv"),
  paste0("- results/", gse, "_pathway_score_tests.csv"),
  paste0("- results/", gse, "_immune_marker_expression.csv"),
  paste0("- figures/Supplementary_Figure_", gse, "_independent_validation.png")
), file.path("results", paste0(gse, "_analysis_summary.md")))

write_session_info()
message("Done: ", gse, " independent validation")
