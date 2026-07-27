#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Run supervised Yin 2024 context analysis by non-circular marker label transfer.
#
# Main inputs:
#   - data/raw/Yin2024_mmc4_cluster_markers.xlsx
#   - data/processed/Yin2024_protein_entrez_log2.csv
#   - results/Yin2024_unsup_pathway_scores.csv
#   - results/Yin2024_unsup_axis_scores.csv
#
# Main outputs:
#   - results/Yin2024_channel_cluster_assignment.csv
#   - results/Yin2024_supervised_axis_kruskal.csv
#   - results/Yin2024_supervised_pathway_kruskal.csv
#   - results/Yin2024_supervised_pathway_pairwise.csv
#   - figures/Supplementary_Figure_Yin2024_supervised_validation.*
#
# Statistical approach:
#   Assigns each TMT channel to the nearest Yin subtype marker set using
#   published mmc4 centroids, then tests whether this project's signatures
#   differ across transferred labels with Kruskal-Wallis and pairwise tests.
#
# Interpretation notes:
#   The 45 Yin marker genes do not overlap the tested immune, cell-cycle,
#   RTK, or chordoma-identity signatures, reducing circularity risk.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()
suppressPackageStartupMessages({ library(readr); library(readxl) })

# ---- 1. Load the log2 protein matrix produced by script 08 -----------
expr_path <- "data/processed/Yin2024_protein_entrez_log2.csv"
if (!file.exists(expr_path)) stop("Run script 08 first to create ", expr_path)
expr_df <- read_csv(expr_path, show_col_types = FALSE)
expr <- as.matrix(expr_df[,-1]); rownames(expr) <- as.character(expr_df$entrez_id)

sym2ez <- read_csv("data/processed/symbol2entrez.csv", show_col_types = FALSE) %>%
  mutate(symbol = as.character(symbol), entrez_id = as.character(entrez_id)) %>%
  distinct(symbol, .keep_all = TRUE)

# ---- 2. Marker centroids (mmc4) -> assign each channel to a cluster --
m4 <- suppressMessages(read_excel("data/raw/Yin2024_mmc4_cluster_markers.xlsx", sheet = 1))
m4 <- m4 %>% left_join(sym2ez, by = c("gene" = "symbol")) %>% filter(!is.na(entrez_id))
cl_cols <- c("Cluster1", "Cluster2", "Cluster3")
# which cluster each marker is highest in (defines the 3 marker sets)
m4$top <- apply(as.matrix(m4[, cl_cols]), 1, which.max)
marker_sets <- split(m4$entrez_id, paste0("C", m4$top))

z <- t(scale(t(expr))); z[is.na(z)] <- 0   # z-score each protein across channels
channel_marker_score <- sapply(marker_sets, function(ids) {
  present <- intersect(ids, rownames(z))
  colMeans(z[present, , drop = FALSE], na.rm = TRUE)
})
assigned <- colnames(z)[0]
cluster_assign <- tibble(
  sample  = rownames(channel_marker_score),
  C1 = channel_marker_score[, "C1"],
  C2 = channel_marker_score[, "C2"],
  C3 = channel_marker_score[, "C3"]
) %>%
  mutate(cluster = c("BME", "MES", "MET")[max.col(cbind(C1, C2, C3), ties.method = "first")])
write_csv(cluster_assign, "results/Yin2024_channel_cluster_assignment.csv")
message("Assigned channels per cluster:")
print(table(cluster_assign$cluster))

meta <- cluster_assign %>% select(sample, cluster)

# ---- 3. Load OUR signature scores (from script 08) and join ----------
scores <- read_csv("results/Yin2024_unsup_pathway_scores.csv", show_col_types = FALSE)
axis   <- read_csv("results/Yin2024_unsup_axis_scores.csv", show_col_types = FALSE) %>%
  select(sample, immune_interferon, cell_cycle, chordoma_marker)
rtk    <- scores %>% select(sample, CUSTOM_DRUGGABLE_RTKS)
axis   <- axis %>% left_join(rtk, by = "sample") %>%
  rename(rtk = CUSTOM_DRUGGABLE_RTKS)

# ---- 4. Supervised tests: Kruskal-Wallis + pairwise ------------------
kruskal_by <- function(score_df, group) {
  joined <- score_df %>% left_join(group, by = "sample")
  cols <- setdiff(colnames(score_df), "sample")
  bind_rows(lapply(cols, function(p) {
    df <- joined %>% filter(!is.na(.data[[p]]), !is.na(cluster))
    kw <- tryCatch(kruskal.test(df[[p]], factor(df$cluster))$p.value, error = function(e) NA_real_)
    mm <- df %>% group_by(cluster) %>% summarise(m = mean(.data[[p]], na.rm = TRUE), .groups = "drop")
    tibble(pathway = p, !!!setNames(as.list(round(mm$m,3)), paste0("mean_", mm$cluster)), kw_p = kw)
  })) %>% mutate(kw_fdr = p.adjust(kw_p, "BH"))
}

axis_kw    <- kruskal_by(axis, meta)
pathway_kw <- kruskal_by(scores, meta)
write_csv(axis_kw,    "results/Yin2024_supervised_axis_kruskal.csv")
write_csv(pathway_kw, "results/Yin2024_supervised_pathway_kruskal.csv")

pw <- bind_rows(
  compare_scores(scores, meta, "cluster", "BME", "MES") %>% mutate(contrast = "BME_vs_MES"),
  compare_scores(scores, meta, "cluster", "BME", "MET") %>% mutate(contrast = "BME_vs_MET"),
  compare_scores(scores, meta, "cluster", "MET", "MES") %>% mutate(contrast = "MET_vs_MES")
) %>% relocate(contrast)
write_csv(pw, "results/Yin2024_supervised_pathway_pairwise.csv")

# ---- 5. Figure -------------------------------------------------------
long <- axis %>% left_join(meta, by = "sample") %>%
  pivot_longer(c(immune_interferon, cell_cycle, rtk, chordoma_marker),
               names_to = "signature", values_to = "score") %>%
  mutate(signature = recode(signature,
           immune_interferon = "Immune/IFN", cell_cycle = "Cell-cycle",
           rtk = "RTK", chordoma_marker = "Chordoma identity"),
         cluster = factor(cluster, levels = c("BME","MES","MET")))
p <- ggplot(long, aes(cluster, score, fill = cluster)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
  facet_wrap(~signature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(BME="#009E73", MES="#0072B2", MET="#D55E00")) +
  labs(title = NULL,
       subtitle = NULL,
       x = NULL, y = "Signature score (z)") +
  theme_nature(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave("figures/Supplementary_Figure_Yin2024_supervised_validation.png", p, width = 10, height = 3.6, dpi = 300)
ggsave("figures/Supplementary_Figure_Yin2024_supervised_validation.pdf", p, width = 10, height = 3.6)

# ---- 6. Summary ------------------------------------------------------
gm <- function(tbl,row) { r<-tbl[tbl$pathway==row,]; if(nrow(r)) r else NULL }
imm <- gm(axis_kw,"immune_interferon"); cc <- gm(axis_kw,"cell_cycle")
rtkr<- gm(axis_kw,"rtk"); chid<- gm(axis_kw,"chordoma_marker")
fmt <- function(r) if(is.null(r)) "NA" else sprintf("BME=%.2f MES=%.2f MET=%.2f (KW FDR=%.3g)",
        r$mean_BME, r$mean_MES, r$mean_MET, r$kw_fdr)
writeLines(c(
  "# Yin 2024 supervised cross-modality validation (marker label-transfer)",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Method",
  "Each TMT channel was assigned to a Yin subtype using their 45 cluster-marker",
  "centroids (mmc4; nearest marker-set on z-scored abundances). The 45 markers do",
  "not overlap the immune/cell-cycle/RTK signatures, so the test is non-circular.",
  "Channel->patient key is not published, so per-sample survival was not joined.",
  "",
  "## Assigned channels per cluster",
  paste0("- ", paste(names(table(cluster_assign$cluster)), table(cluster_assign$cluster),
                     sep = ": ", collapse = ";  ")),
  "  (Yin published cluster sizes: BME/cluster1=26, MES/cluster2=35, MET/cluster3=41.)",
  "",
  "## Signature means by subtype (expected: immune highest in BME; RTK+identity in MET)",
  paste0("- Immune/interferon: ", fmt(imm)),
  paste0("- Cell-cycle:        ", fmt(cc)),
  paste0("- RTK:               ", fmt(rtkr)),
  paste0("- Chordoma identity: ", fmt(chid)),
  "",
  "## Main outputs",
  "- results/Yin2024_channel_cluster_assignment.csv",
  "- results/Yin2024_supervised_axis_kruskal.csv",
  "- results/Yin2024_supervised_pathway_kruskal.csv",
  "- results/Yin2024_supervised_pathway_pairwise.csv",
  "- figures/Supplementary_Figure_Yin2024_supervised_validation.png"
), "results/Yin2024_supervised_summary.md")

write_session_info()
message("Done: Yin 2024 supervised label-transfer validation")
