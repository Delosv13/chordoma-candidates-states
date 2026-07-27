#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Run unsupervised cross-modality context analysis in Yin 2024 proteomic channels.
#
# Main inputs:
#   - data/raw/Yin2024_protein_abundance.csv
#   - data/processed/symbol2entrez_primary.csv
#   - data/raw/h.all.v2025.1.Hs.entrez.gmt
#
# Main outputs:
#   - data/processed/Yin2024_protein_entrez_log2.csv
#   - results/Yin2024_unsup_pathway_scores.csv
#   - results/Yin2024_unsup_axis_scores.csv
#   - results/Yin2024_unsup_immune_marker_correlations.csv
#   - results/Yin2024_unsup_cluster_profiles.csv
#   - figures/Supplementary_Figure_Yin2024_unsupervised_validation.*
#
# Statistical approach:
#   Excludes reference-pool channels, maps protein symbols to Entrez, log2
#   transforms protein abundance, scores the same signatures as RNA cohorts,
#   uses Spearman correlations, and performs unsupervised k=3 clustering.
#
# Interpretation notes:
#   Protein columns are TMT channels rather than patient IDs. No published
#   subtype labels are attached in this script.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()
suppressPackageStartupMessages(library(readr))

prot_path <- "data/raw/Yin2024_protein_abundance.csv"
prot <- read_csv(prot_path, show_col_types = FALSE)
gene_col <- colnames(prot)[1]

# Drop reference-pool channels; keep sample TMT channels only
sample_cols <- setdiff(colnames(prot), gene_col)
sample_cols <- sample_cols[!grepl("RefInt", sample_cols, ignore.case = TRUE)]
message("Sample TMT channels: ", length(sample_cols))

# Row IDs are "UniProt_GENESYMBOL"; take symbol after the first underscore
prot$symbol <- sub("^[^_]*_", "", as.character(prot[[gene_col]]))

# symbol -> Entrez map built from NCBI Homo_sapiens.gene_info (primary + synonyms)
sym2ez <- read_csv("data/processed/symbol2entrez.csv", show_col_types = FALSE) %>%
  mutate(symbol = as.character(symbol), entrez_id = as.character(entrez_id)) %>%
  filter(!is.na(entrez_id), entrez_id != "") %>%
  distinct(symbol, .keep_all = TRUE)

mapped <- prot %>%
  inner_join(sym2ez, by = "symbol") %>%
  select(entrez_id, all_of(sample_cols)) %>%
  group_by(entrez_id) %>%
  summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
expr <- matrix_from_gene_table(mapped, "entrez_id", sample_cols)
message("Mapped proteins to Entrez: ", nrow(expr))

# Raw intensities -> log2
expr <- log2(pmax(expr, 0) + 1)
write_csv(as.data.frame(expr) %>% rownames_to_column("entrez_id"),
          "data/processed/Yin2024_protein_entrez_log2.csv")

# ---- Score the SAME gene sets / axes ---------------------------------
gene_sets <- selected_gene_sets()
scores <- score_gene_sets(expr, gene_sets)
axis_scores <- axis_scores_from_pathways(scores)
write_csv(scores, "results/Yin2024_unsup_pathway_scores.csv")
write_csv(axis_scores, "results/Yin2024_unsup_axis_scores.csv")
write_csv(gene_set_overlap(expr, gene_sets), "results/Yin2024_unsup_pathway_gene_overlap.csv")

# ---- (1) immune-marker proteins vs immune axis -----------------------
immune_markers <- tibble(
  gene_symbol = c("CD274", "PDCD1", "HAVCR2", "STAT1", "CXCL9", "GZMB", "CD8A",
                  "CD3G", "CD14", "CXCL13"),
  entrez_id   = c("29126", "5133", "84868", "6772", "4283", "3002", "925",
                  "917", "929", "10563")
)
imm_axis <- setNames(axis_scores$immune_interferon, axis_scores$sample)
marker_cor <- lapply(seq_len(nrow(immune_markers)), function(i) {
  gid <- immune_markers$entrez_id[i]
  if (!gid %in% rownames(expr)) return(tibble(gene_symbol = immune_markers$gene_symbol[i],
                                               entrez_id = gid, n_obs = 0L, rho = NA_real_))
  v <- as.numeric(expr[gid, axis_scores$sample])
  tibble(gene_symbol = immune_markers$gene_symbol[i], entrez_id = gid,
         n_obs = sum(is.finite(v)),
         rho = if (sum(is.finite(v)) >= 5)
           suppressWarnings(cor(v, imm_axis, method = "spearman", use = "complete.obs")) else NA_real_)
}) %>% bind_rows()
write_csv(marker_cor, "results/Yin2024_unsup_immune_marker_correlations.csv")

# ---- (2) immune<->cell-cycle relationship (vs RNA cohorts) -----------
rho_yin <- cor(axis_scores$immune_interferon, axis_scores$cell_cycle, method = "spearman")
g <- read_csv("results/GSE230168_axis_scores.csv", show_col_types = FALSE)
gt <- g[g$type != "nucleus_pulposus", ]
rho_230168 <- cor(gt$immune_interferon, gt$cell_cycle, method = "spearman")

# ---- (3) unsupervised k=3 clustering on the scored signatures --------
feat <- scores %>% column_to_rownames("sample") %>% as.matrix()
feat[is.na(feat)] <- 0
set.seed(123)
km <- kmeans(scale(feat), centers = 3, nstart = 50)
axis_scores$cluster <- factor(paste0("k", km$cluster))
clust_profile <- axis_scores %>%
  group_by(cluster) %>%
  summarise(n = n(),
            immune = mean(immune_interferon, na.rm = TRUE),
            cell_cycle = mean(cell_cycle, na.rm = TRUE),
            chordoma_marker = mean(chordoma_marker, na.rm = TRUE),
            .groups = "drop")
# add RTK mean per cluster
rtk <- scores %>% select(sample, CUSTOM_DRUGGABLE_RTKS) %>%
  left_join(axis_scores %>% select(sample, cluster), by = "sample") %>%
  group_by(cluster) %>% summarise(rtk = mean(CUSTOM_DRUGGABLE_RTKS, na.rm = TRUE), .groups = "drop")
clust_profile <- clust_profile %>% left_join(rtk, by = "cluster")
write_csv(axis_scores, "results/Yin2024_unsup_axis_scores_clustered.csv")
write_csv(clust_profile, "results/Yin2024_unsup_cluster_profiles.csv")

# ---- Figure ----------------------------------------------------------
p <- ggplot(axis_scores, aes(immune_interferon, cell_cycle, color = cluster)) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_color_manual(values = unname(npg_palette())) +
  labs(title = NULL,
       subtitle = NULL,
       x = "Immune/interferon score (protein)", y = "Cell-cycle score (protein)") +
  theme_nature(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("figures/Supplementary_Figure_Yin2024_unsupervised_validation.png", p, width = 6.8, height = 5, dpi = 300)
ggsave("figures/Supplementary_Figure_Yin2024_unsupervised_validation.pdf", p, width = 6.8, height = 5)

# ---- Summary ---------------------------------------------------------
writeLines(c(
  "# Yin 2024 cross-modality validation (unsupervised)",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## What this tests (no subtype labels needed)",
  "Reprojection of the RNA-defined state signatures onto the Yin 2024 protein",
  "matrix (iProX IPX0009051000), testing internal consistency and structure.",
  "",
  "## Key numbers",
  paste0("- Sample TMT channels scored: ", nrow(axis_scores), " (reference pools dropped)."),
  paste0("- Proteins mapped to Entrez: ", nrow(expr), "."),
  paste0("- immune~cell-cycle Spearman rho (protein): ", round(rho_yin, 3),
         "  vs GSE230168 tumors (RNA): ", round(rho_230168, 3), "."),
  "- Immune-marker protein correlations with the immune axis: see",
  "  results/Yin2024_unsup_immune_marker_correlations.csv",
  "- Unsupervised k=3 cluster profiles: results/Yin2024_unsup_cluster_profiles.csv",
  "",
  "## Interpretation",
  "- Immune-marker proteins (CD8A/CD274/STAT1/CXCL9/CXCL13...) tracking the",
  "  immune axis = the immune state exists at the protein level.",
  "- One k=3 cluster high on immune + one high on chordoma_marker/RTK would",
  "  recapitulate Yin's BM-dominant vs MET-mediated subtypes from our signatures.",
  "- The formal supervised directional test is implemented in script 09 via marker label-transfer.",
  "",
  "## Main outputs",
  "- results/Yin2024_unsup_immune_marker_correlations.csv",
  "- results/Yin2024_unsup_cluster_profiles.csv",
  "- results/Yin2024_unsup_axis_scores_clustered.csv",
  "- figures/Supplementary_Figure_Yin2024_unsupervised_validation.png"
), "results/Yin2024_unsupervised_summary.md")

write_session_info()
message("Done: Yin 2024 unsupervised cross-modality validation")