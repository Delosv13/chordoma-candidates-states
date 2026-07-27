#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Run the primary held-out naive RNA validation in GSE312699.
#
# Main inputs:
#   - data/raw/GSE312699/GSE312699_series_matrix.txt.gz
#   - data/raw/GSE312699/GPL15207_family.soft.gz
#   - data/raw/h.all.v2025.1.Hs.entrez.gmt
#
# Main outputs:
#   - data/processed/GSE312699_log2expr_entrez.csv
#   - data/processed/GSE312699_metadata_curated.csv
#   - results/GSE312699_pathway_scores.csv
#   - results/GSE312699_axis_scores.csv
#   - results/GSE312699_pathway_score_tests.csv
#   - results/GSE312699_axis_score_tests.csv
#   - results/GSE312699_primary_vs_recurrent_signature.csv
#   - figures/Supplementary_Figure_GSE312699_naive_rna_validation.*
#
# Statistical approach:
#   GEO Series Matrix values are used as normalized microarray expression.
#   GPL15207 probe annotations are collapsed to Entrez genes by mean probe
#   expression. Prespecified signatures are scored as mean sample-wise z-scores.
#   Primary-vs-recurrent comparisons use paired Wilcoxon tests by patient pair,
#   and the gene-level signature uses limma with patient pair as a blocking
#   covariate.
#
# Interpretation notes:
#   GSE312699 is the primary held-out naive RNA validation cohort in the revised
#   evidence hierarchy. The analysis is cross-platform and small-n (six pairs),
#   so interpretation should emphasize direction and confidence boundaries.
# -------------------------------------------------------------------------

source("scripts/00_chordoma_utils.R")
ensure_project_dirs()

raw_dir <- "data/raw/GSE312699"
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

series_path <- file.path(raw_dir, "GSE312699_series_matrix.txt.gz")
platform_path <- file.path(raw_dir, "GPL15207_family.soft.gz")
platform_cache <- "data/processed/GSE312699_GPL15207_probe_entrez_map.csv"
series_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312699/matrix/GSE312699_series_matrix.txt.gz"
platform_url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15207/soft/GPL15207_family.soft.gz"

download_if_missing <- function(url, dest) {
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    message("Downloading ", basename(dest), " ...")
    utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  }
}

download_if_missing(series_url, series_path)
download_if_missing(platform_url, platform_path)

strip_quotes <- function(x) {
  gsub('^"|"$', "", x)
}

split_quoted_tsv <- function(line) {
  strip_quotes(strsplit(line, "\t", fixed = TRUE)[[1]])
}

extract_series_matrix <- function(path) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)

  sample_line <- function(prefix) {
    hit <- grep(paste0("^", prefix, "\t"), lines, value = TRUE)
    if (length(hit) == 0) {
      stop("Missing Series Matrix metadata line: ", prefix)
    }
    split_quoted_tsv(hit[[1]])[-1]
  }

  titles <- sample_line("!Sample_title")
  accessions <- sample_line("!Sample_geo_accession")
  characteristics <- sample_line("!Sample_characteristics_ch1")

  meta <- tibble(
    sample = accessions,
    title = titles,
    characteristics = characteristics
  ) %>%
    mutate(
      patient_pair = sub("^.*chordoma_([0-9]+)_(primary|recurrence).*$", "\\1", title),
      condition = if_else(grepl("recurrence", title, ignore.case = TRUE), "recurrent", "primary"),
      condition = factor(condition, levels = c("primary", "recurrent")),
      patient_pair = factor(patient_pair)
    )

  begin <- grep("^!series_matrix_table_begin", lines)
  end <- grep("^!series_matrix_table_end", lines)
  if (length(begin) != 1 || length(end) != 1 || end <= begin) {
    stop("Could not locate the GSE312699 expression table in the Series Matrix.")
  }

  table_lines <- lines[(begin + 1):(end - 1)]
  expr_df <- read_tsv(
    I(paste(table_lines, collapse = "\n")),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    rename(probe_id = ID_REF)

  list(meta = meta, expr_df = expr_df)
}

extract_platform_annotation <- function(path) {
  lines <- readLines(gzfile(path, "rt"), warn = FALSE)
  begin <- grep("^!platform_table_begin", lines)
  if (length(begin) != 1) {
    stop("Could not locate GPL15207 platform table.")
  }

  after <- lines[(begin + 1):length(lines)]
  stop_at <- which(startsWith(after, "!"))
  if (length(stop_at) > 0) {
    after <- after[seq_len(stop_at[[1]] - 1)]
  }
  after <- after[nzchar(after)]

  annot <- read_tsv(
    I(paste(after, collapse = "\n")),
    show_col_types = FALSE,
    progress = FALSE
  )

  required <- c("ID", "Gene Symbol", "Entrez Gene")
  missing <- setdiff(required, colnames(annot))
  if (length(missing) > 0) {
    stop("Missing GPL15207 annotation column(s): ", paste(missing, collapse = ", "))
  }

  annot %>%
    transmute(
      probe_id = as.character(ID),
      gene_symbol = as.character(`Gene Symbol`),
      entrez_raw = as.character(`Entrez Gene`)
    ) %>%
    filter(!is.na(entrez_raw), entrez_raw != "", entrez_raw != "---") %>%
    separate_rows(entrez_raw, sep = "\\s*///\\s*|\\s*//\\s*|\\s*;\\s*|\\s*,\\s*") %>%
    mutate(entrez_id = gsub("[^0-9]", "", entrez_raw)) %>%
    filter(entrez_id != "") %>%
    distinct(probe_id, entrez_id, .keep_all = TRUE)
}

paired_score_test <- function(scores, meta, group_col = "condition",
                              case_value = "recurrent", control_value = "primary",
                              pair_col = "patient_pair") {
  joined <- scores %>% left_join(meta, by = "sample")
  score_cols <- setdiff(colnames(scores), "sample")

  bind_rows(lapply(score_cols, function(pathway) {
    wide <- joined %>%
      select(all_of(c(pair_col, group_col, pathway))) %>%
      filter(.data[[group_col]] %in% c(case_value, control_value)) %>%
      pivot_wider(names_from = all_of(group_col), values_from = all_of(pathway)) %>%
      filter(is.finite(.data[[case_value]]), is.finite(.data[[control_value]]))

    delta <- wide[[case_value]] - wide[[control_value]]
    p <- if (length(delta) >= 2 && sd(delta, na.rm = TRUE) > 0) {
      suppressWarnings(wilcox.test(wide[[case_value]], wide[[control_value]],
                                   paired = TRUE, exact = FALSE)$p.value)
    } else {
      NA_real_
    }

    tibble(
      pathway = pathway,
      group_a = case_value,
      group_b = control_value,
      n_pairs = nrow(wide),
      mean_a = mean(wide[[case_value]], na.rm = TRUE),
      mean_b = mean(wide[[control_value]], na.rm = TRUE),
      delta_a_minus_b = mean(delta, na.rm = TRUE),
      median_pair_delta = median(delta, na.rm = TRUE),
      p_value = p
    )
  })) %>%
    mutate(fdr = p.adjust(p_value, method = "BH"))
}

message("Reading GSE312699 Series Matrix")
series <- extract_series_matrix(series_path)
meta <- series$meta
expr_df <- series$expr_df
sample_cols <- meta$sample

if (file.exists(platform_cache)) {
  message("Reading cached GPL15207 probe-to-Entrez map")
  probe_map <- read_csv(platform_cache, show_col_types = FALSE)
} else {
  message("Reading GPL15207 annotation")
  probe_map <- extract_platform_annotation(platform_path)
  write_csv(probe_map, platform_cache)
}

expr_mapped <- expr_df %>%
  inner_join(probe_map, by = "probe_id") %>%
  select(entrez_id, all_of(sample_cols)) %>%
  group_by(entrez_id) %>%
  summarise(across(all_of(sample_cols), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

log_expr <- matrix_from_gene_table(expr_mapped, "entrez_id", sample_cols)
log_expr <- log_expr[rowSums(is.finite(log_expr)) == ncol(log_expr), , drop = FALSE]
log_expr <- log_expr[apply(log_expr, 1, stats::var, na.rm = TRUE) > 0, , drop = FALSE]

gene_sets <- selected_gene_sets()
scores <- score_gene_sets(log_expr, gene_sets)
axis_scores <- axis_scores_from_pathways(scores)

pathway_tests <- paired_score_test(scores, meta)
axis_tests <- paired_score_test(axis_scores, meta)
signature <- differential_signature_limma(log_expr, meta, "condition", "recurrent", "primary",
                                          block_col = "patient_pair")

write_csv(meta, "data/processed/GSE312699_metadata_curated.csv")
write_csv(as.data.frame(log_expr) %>% rownames_to_column("entrez_id"),
          "data/processed/GSE312699_log2expr_entrez.csv")
write_csv(scores %>% left_join(meta, by = "sample"), "results/GSE312699_pathway_scores.csv")
write_csv(axis_scores %>% left_join(meta, by = "sample"), "results/GSE312699_axis_scores.csv")
write_csv(gene_set_overlap(log_expr, gene_sets), "results/GSE312699_pathway_gene_overlap.csv")
write_csv(pathway_tests, "results/GSE312699_pathway_score_tests.csv")
write_csv(axis_tests, "results/GSE312699_axis_score_tests.csv")
write_csv(signature, "results/GSE312699_primary_vs_recurrent_signature.csv")

labels <- pathway_labels()
selected <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "CUSTOM_ANTIGEN_PRESENTATION",
  "CUSTOM_CHORDOMA_MARKERS",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "CUSTOM_CDK46_AXIS",
  "CUSTOM_PLK1_AXIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS"
)

delta_df <- pathway_tests %>%
  filter(pathway %in% selected) %>%
  mutate(pathway_label = factor(labels[pathway], levels = rev(labels[selected])))

p_delta <- ggplot(delta_df, aes(x = delta_a_minus_b, y = pathway_label)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  geom_col(fill = "#0072B2", color = "grey20", linewidth = 0.2, width = 0.72) +
  labs(
    title = "a",
    x = "Mean paired score difference (recurrent - primary)",
    y = NULL
  ) +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

axis_long <- axis_scores %>%
  left_join(meta, by = "sample") %>%
  pivot_longer(c(immune_interferon, cell_cycle, chordoma_marker),
               names_to = "axis", values_to = "score") %>%
  mutate(axis = recode(axis,
                       immune_interferon = "Immune/interferon",
                       cell_cycle = "Cell-cycle",
                       chordoma_marker = "Chordoma markers"))

p_axis <- ggplot(axis_long, aes(x = condition, y = score, group = patient_pair)) +
  geom_line(color = "grey55", linewidth = 0.35, alpha = 0.85) +
  geom_point(aes(fill = condition), shape = 21, size = 2.4, color = "grey20") +
  facet_wrap(~axis, scales = "free_y") +
  scale_fill_manual(values = c(primary = "#0072B2", recurrent = "#D55E00")) +
  labs(title = "b", x = NULL, y = "Mean z-score") +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

plot_genes <- c(TBXT = "6862", CD274 = "29126", HLA.B = "3106", TAP1 = "6890",
                MKI67 = "4288", CDK4 = "1019", PLK1 = "5347", ASNS = "440")
gene_df <- bind_rows(lapply(names(plot_genes), function(symbol) {
  gid <- plot_genes[[symbol]]
  if (!gid %in% rownames(log_expr)) return(tibble())
  tibble(sample = colnames(log_expr), gene = gsub("\\.", "-", symbol),
         entrez_id = gid, expression = as.numeric(log_expr[gid, ]))
})) %>%
  left_join(meta, by = "sample")

p_gene <- ggplot(gene_df, aes(x = condition, y = expression, group = patient_pair)) +
  geom_line(color = "grey60", linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(fill = condition), shape = 21, size = 2.1, color = "grey20") +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(primary = "#0072B2", recurrent = "#D55E00")) +
  labs(title = "c", x = NULL, y = "Normalized expression") +
  theme_nature(base_size = 8.8) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

volcano_df <- signature %>%
  mutate(
    neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
    label = if_else(entrez_id %in% plot_genes,
                    names(plot_genes)[match(entrez_id, plot_genes)],
                    NA_character_)
  )

p_volcano <- ggplot(volcano_df, aes(x = log2fc, y = neg_log10_fdr)) +
  geom_point(color = "grey65", size = 0.75, alpha = 0.45) +
  geom_point(data = filter(volcano_df, !is.na(label)), color = "#D55E00", size = 2) +
  ggrepel::geom_text_repel(aes(label = gsub("\\.", "-", label)), max.overlaps = 20,
                           size = 3, na.rm = TRUE) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  labs(title = "d",
       x = "limma log2 fold-change", y = "-log10(FDR)") +
  theme_nature(base_size = 10) +
  theme(panel.grid.minor = element_blank())

fig <- plot_grid(
  plot_grid(p_delta, p_axis, nrow = 1, rel_widths = c(1.05, 1)),
  plot_grid(p_gene, p_volcano, nrow = 1, rel_widths = c(1.25, 1)),
  ncol = 1
)

ggsave("figures/Supplementary_Figure_GSE312699_naive_rna_validation.png",
       fig, width = 13, height = 8.8, dpi = 300)
ggsave("figures/Supplementary_Figure_GSE312699_naive_rna_validation.pdf",
       fig, width = 13, height = 8.8)

axis_brief <- axis_tests %>%
  mutate(line = paste0(
    "- ", pathway, ": delta = ", round(delta_a_minus_b, 3),
    ", paired Wilcoxon p = ", signif(p_value, 3),
    ", FDR = ", signif(fdr, 3)
  )) %>%
  pull(line)

writeLines(c(
  "# GSE312699 primary held-out naive RNA validation",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "## Scope",
  "Primary-vs-recurrent paired microarray validation in six matched human chordoma patients.",
  "",
  "## Basic checks",
  paste0("- Samples: ", nrow(meta), " (", n_distinct(meta$patient_pair), " matched pairs)."),
  paste0("- Series Matrix probes: ", nrow(expr_df), "."),
  paste0("- GPL15207 annotated probe-to-Entrez rows used: ", nrow(probe_map), "."),
  paste0("- Collapsed Entrez genes retained: ", nrow(log_expr), "."),
  "",
  "## Axis-level paired tests",
  axis_brief,
  "",
  "## Interpretation guide",
  "- This is the primary held-out naive RNA validation cohort in the revised evidence hierarchy.",
  "- Because n = 6 matched pairs and the platform is microarray, use this as directional validation of prespecified axes rather than a discovery source.",
  "- PXD025894/Shen 2021 remains a supporting targeted proteomic AvsB check, not a proteome-wide primary validation.",
  "",
  "## Main outputs",
  "- data/processed/GSE312699_GPL15207_probe_entrez_map.csv",
  "- data/processed/GSE312699_log2expr_entrez.csv",
  "- data/processed/GSE312699_metadata_curated.csv",
  "- results/GSE312699_pathway_scores.csv",
  "- results/GSE312699_axis_scores.csv",
  "- results/GSE312699_pathway_score_tests.csv",
  "- results/GSE312699_axis_score_tests.csv",
  "- results/GSE312699_primary_vs_recurrent_signature.csv",
  "- figures/Supplementary_Figure_GSE312699_naive_rna_validation.png"
), "results/GSE312699_analysis_summary.md")

write_session_info()
message("Done: GSE312699 naive RNA validation")
