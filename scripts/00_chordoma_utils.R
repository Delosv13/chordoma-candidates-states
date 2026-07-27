# -------------------------------------------------------------------------
# Methodological header
# Purpose:
#   Shared utility functions for the chordoma reproducibility pipeline.
#
# Main inputs:
#   - Public expression matrices and metadata loaded by downstream scripts.
#   - MSigDB Hallmark GMT and NCBI/HGNC mapping resources.
#
# Main outputs:
#   - No direct analysis outputs when sourced alone.
#   - Downstream scripts use these functions to write processed matrices,
#     score tables, statistical summaries, figures, and session metadata.
#
# Statistical approach:
#   Provides Entrez harmonization, log2(CPM+1) transformation, mean
#   sample-wise z-score pathway scoring, Wilcoxon/BH-FDR group comparisons,
#   legacy Welch t-test signatures, and canonical limma-trend/eBayes
#   differential signatures for small-n log-expression contrasts.
#
# Interpretation notes:
#   This file centralizes reusable methods so the statistical definitions are
#   not duplicated across datasets. Cross-study inference remains at the
#   signature/pathway and direction-of-effect level; matrices are not
#   batch-merged here.
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(cowplot)
  library(limma)
})

ensure_project_dirs <- function() {
  dirs <- c("data/raw", "data/processed", "results", "figures", "scripts")
  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parsed <- strsplit(lines, "\t", fixed = TRUE)
  gene_sets <- lapply(parsed, function(x) unique(as.character(x[-c(1, 2)])))
  names(gene_sets) <- vapply(parsed, function(x) x[[1]], character(1))
  gene_sets
}

read_human_gene2ensembl <- function(path = "data/raw/NCBI_gene2ensembl.gz",
                                    cache_path = "data/processed/NCBI_human_gene2ensembl.csv",
                                    chunk_size = 200000) {
  if (file.exists(cache_path)) {
    return(read_csv(cache_path, show_col_types = FALSE))
  }

  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- sub("^#", "", readLines(con, n = 1, warn = FALSE))
  cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  chunks <- list()
  i <- 1L

  repeat {
    lines <- readLines(con, n = chunk_size, warn = FALSE)
    if (length(lines) == 0) break
    lines <- lines[startsWith(lines, "9606\t")]
    if (length(lines) == 0) next
    chunks[[i]] <- read.table(
      text = paste(lines, collapse = "\n"),
      sep = "\t",
      header = FALSE,
      col.names = cols,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    )
    i <- i + 1L
  }

  if (length(chunks) == 0) {
    stop("No human tax_id=9606 records found in NCBI gene2ensembl file.")
  }

  gene_map <- bind_rows(chunks) %>%
    transmute(
      entrez_id = as.character(GeneID),
      ensembl_gene_id = sub("\\..*$", "", Ensembl_gene_identifier)
    ) %>%
    filter(
      !is.na(ensembl_gene_id),
      !is.na(entrez_id),
      ensembl_gene_id != "-",
      entrez_id != "-"
    ) %>%
    distinct(ensembl_gene_id, entrez_id)

  write_csv(gene_map, cache_path)
  gene_map
}

map_ensembl_matrix_to_entrez <- function(df, gene_col, sample_cols, gene_map) {
  df %>%
    mutate(ensembl_gene_id = sub("\\..*$", "", .data[[gene_col]])) %>%
    inner_join(gene_map, by = "ensembl_gene_id") %>%
    select(entrez_id, all_of(sample_cols)) %>%
    group_by(entrez_id) %>%
    summarise(across(everything(), ~sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
}

matrix_from_gene_table <- function(df, gene_col, sample_cols) {
  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  rownames(mat) <- as.character(df[[gene_col]])
  storage.mode(mat) <- "numeric"
  mat
}

log2_cpm <- function(counts, min_cpm = 1, min_samples = 3) {
  counts[is.na(counts)] <- 0
  lib_size <- colSums(counts, na.rm = TRUE)
  cpm <- sweep(counts, 2, lib_size / 1e6, "/")
  keep <- rowSums(cpm > min_cpm) >= min_samples
  log2(cpm[keep, , drop = FALSE] + 1)
}

collapse_entrez_matrix <- function(expr) {
  df <- as.data.frame(expr) %>% rownames_to_column("entrez_id")
  df %>%
    group_by(entrez_id) %>%
    summarise(across(everything(), ~mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("entrez_id") %>%
    as.matrix()
}

selected_gene_sets <- function(gmt_path = "data/raw/h.all.v2025.1.Hs.entrez.gmt") {
  hallmark <- read_gmt(gmt_path)
  selected_sets <- c(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
  )

  custom_sets <- list(
    CUSTOM_ANTIGEN_PRESENTATION = c("567", "3105", "3106", "3107", "3122", "3123", "3113", "3115", "972", "6890", "6891", "5696", "5698", "5699", "4261", "84166"),
    CUSTOM_CHORDOMA_MARKERS = c("6862", "3856", "3875", "3880", "4072", "100133941", "6662"),
    CUSTOM_CDK46_AXIS = c("1019", "1021", "595", "894", "896", "5925", "1869", "1870", "1871", "1029", "1030"),
    CUSTOM_PLK1_AXIS = c("5347", "11065", "991", "9133", "1111", "3832", "4751", "898", "4171", "4172"),
    CUSTOM_DRUGGABLE_RTKS = c("1956", "5156", "5159", "3791", "3815", "4233", "2260", "2263", "2261")
  )

  gene_sets <- c(hallmark[selected_sets], custom_sets)
  gene_sets[!vapply(gene_sets, is.null, logical(1))]
}

pathway_labels <- function() {
  c(
    HALLMARK_INTERFERON_ALPHA_RESPONSE = "Interferon alpha",
    HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon gamma",
    HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6-JAK-STAT3",
    HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
    HALLMARK_COMPLEMENT = "Complement",
    HALLMARK_E2F_TARGETS = "E2F targets",
    HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
    HALLMARK_MTORC1_SIGNALING = "mTORC1 signaling",
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
    HALLMARK_ANGIOGENESIS = "Angiogenesis",
    HALLMARK_WNT_BETA_CATENIN_SIGNALING = "WNT/beta-catenin",
    CUSTOM_ANTIGEN_PRESENTATION = "Antigen presentation",
    CUSTOM_CHORDOMA_MARKERS = "Chordoma markers",
    CUSTOM_CDK46_AXIS = "CDK4/6 axis",
    CUSTOM_PLK1_AXIS = "PLK1 axis",
    CUSTOM_DRUGGABLE_RTKS = "Druggable RTKs"
  )
}

score_gene_sets <- function(log_expr, gene_sets) {
  z_expr <- t(scale(t(log_expr)))
  z_expr[is.na(z_expr)] <- 0

  score_one <- function(genes) {
    present <- intersect(as.character(genes), rownames(z_expr))
    if (length(present) < 3) {
      return(rep(NA_real_, ncol(z_expr)))
    }
    colMeans(z_expr[present, , drop = FALSE], na.rm = TRUE)
  }

  score_matrix <- do.call(rbind, lapply(gene_sets, score_one))
  colnames(score_matrix) <- colnames(z_expr)
  as.data.frame(t(score_matrix)) %>% rownames_to_column("sample")
}

gene_set_overlap <- function(log_expr, gene_sets) {
  tibble(
    pathway = names(gene_sets),
    n_genes_requested = vapply(gene_sets, length, integer(1)),
    n_genes_present = vapply(gene_sets, function(x) length(intersect(as.character(x), rownames(log_expr))), integer(1))
  )
}

axis_scores_from_pathways <- function(scores) {
  scores %>%
    transmute(
      sample,
      immune_interferon = rowMeans(across(any_of(c(
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "CUSTOM_ANTIGEN_PRESENTATION"
      ))), na.rm = TRUE),
      cell_cycle = rowMeans(across(any_of(c(
        "HALLMARK_E2F_TARGETS",
        "HALLMARK_G2M_CHECKPOINT",
        "CUSTOM_CDK46_AXIS",
        "CUSTOM_PLK1_AXIS"
      ))), na.rm = TRUE),
      chordoma_marker = CUSTOM_CHORDOMA_MARKERS
    )
}

compare_scores <- function(scores, meta, group_col, group_a, group_b) {
  joined <- scores %>% left_join(meta, by = "sample")
  score_cols <- setdiff(colnames(scores), "sample")
  bind_rows(lapply(score_cols, function(pathway) {
    sub <- joined %>% filter(.data[[group_col]] %in% c(group_a, group_b))
    x <- sub %>% filter(.data[[group_col]] == group_a) %>% pull(all_of(pathway))
    y <- sub %>% filter(.data[[group_col]] == group_b) %>% pull(all_of(pathway))
    p <- if (sum(!is.na(x)) >= 2 && sum(!is.na(y)) >= 2) {
      suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value)
    } else {
      NA_real_
    }
    mean_a <- mean(x, na.rm = TRUE)
    mean_b <- mean(y, na.rm = TRUE)
    tibble(
      pathway = pathway,
      group_a = group_a,
      group_b = group_b,
      n_a = sum(!is.na(x)),
      n_b = sum(!is.na(y)),
      mean_a = mean_a,
      mean_b = mean_b,
      delta_a_minus_b = mean_a - mean_b,
      p_value = p
    )
  })) %>%
    mutate(fdr = p.adjust(p_value, method = "BH"))
}

differential_signature <- function(log_expr, meta, group_col, case_value, control_value) {
  common <- intersect(colnames(log_expr), meta$sample)
  log_expr <- log_expr[, common, drop = FALSE]
  meta <- meta %>% filter(sample %in% common) %>% arrange(match(sample, common))
  case_samples <- meta %>% filter(.data[[group_col]] == case_value) %>% pull(sample)
  control_samples <- meta %>% filter(.data[[group_col]] == control_value) %>% pull(sample)

  out <- lapply(rownames(log_expr), function(g) {
    x <- as.numeric(log_expr[g, case_samples, drop = TRUE])
    y <- as.numeric(log_expr[g, control_samples, drop = TRUE])
    p <- if (sum(!is.na(x)) >= 2 && sum(!is.na(y)) >= 2) {
      suppressWarnings(t.test(x, y)$p.value)
    } else {
      NA_real_
    }
    tibble(
      entrez_id = g,
      mean_case = mean(x, na.rm = TRUE),
      mean_control = mean(y, na.rm = TRUE),
      log2fc = mean_case - mean_control,
      p_value = p
    )
  })

  bind_rows(out) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    arrange(desc(abs(log2fc)))
}

# Versao moderada (limma-trend + empirical Bayes) da assinatura diferencial.
# Preferida sobre o t-test de Welch no regime de n pequeno deste estudo:
# eBayes(trend=TRUE) encolhe a variancia de cada gene em direcao a uma tendencia
# media-variancia comum (Smyth 2004), e block_col permite bloquear o modelo PDX
# (CD3/CD7) em GSE216418 em vez de confundi-lo com o tratamento. Entrada esperada:
# log2(TPM+1)/log2(CPM+1) -> usa-se limma-trend (nao voom). Colunas de saida
# compativeis com differential_signature() (entrez_id, mean_*, log2fc, p_value, fdr).
differential_signature_limma <- function(log_expr, meta, group_col,
                                         case_value, control_value,
                                         block_col = NULL) {
  common <- intersect(colnames(log_expr), meta$sample)
  log_expr <- log_expr[, common, drop = FALSE]
  meta <- meta[match(common, meta$sample), , drop = FALSE]

  grp <- factor(ifelse(meta[[group_col]] == case_value, "case", "control"),
                levels = c("control", "case"))
  if (!is.null(block_col) && block_col %in% names(meta) &&
      dplyr::n_distinct(meta[[block_col]]) > 1) {
    block <- factor(meta[[block_col]])
    design <- model.matrix(~ block + grp)
  } else {
    design <- model.matrix(~ grp)
  }

  fit <- lmFit(log_expr, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
  tt <- tt[rownames(log_expr), , drop = FALSE]

  tibble(
    entrez_id    = rownames(log_expr),
    mean_case    = rowMeans(log_expr[, grp == "case",    drop = FALSE], na.rm = TRUE),
    mean_control = rowMeans(log_expr[, grp == "control", drop = FALSE], na.rm = TRUE),
    log2fc       = tt$logFC,
    avg_expr     = tt$AveExpr,
    mod_t        = tt$t,
    p_value      = tt$P.Value,
    fdr          = tt$adj.P.Val
  ) %>% arrange(desc(abs(log2fc)))
}

signature_set_scores <- function(log_expr, signature, prefix, top_n = 200) {
  sig <- signature %>% filter(!is.na(log2fc))
  up <- sig %>% arrange(desc(log2fc)) %>% slice_head(n = top_n) %>% pull(entrez_id)
  down <- sig %>% arrange(log2fc) %>% slice_head(n = top_n) %>% pull(entrez_id)
  score_gene_sets(log_expr, setNames(list(up, down), paste0(prefix, c("_UP", "_DOWN"))))
}

write_session_info <- function(path = "results/session_info.txt") {
  capture.output(sessionInfo(), file = path)
}

# -------------------------------------------------------------------------
# Color scheme (publication-grade, evidence-based and color-vision-safe)
#
# Rationale, with primary references:
#   Categorical colors use the Okabe-Ito "Color Universal Design" palette,
#   the de-facto standard for colorblind-safe scientific figures, popularized
#   for biology by:
#     Wong, B. (2011) Points of view: Color blindness. Nature Methods 8:441.
#       https://doi.org/10.1038/nmeth.1618
#   These eight hues were verified by Okabe & Ito on color-vision-deficiency
#   simulators to stay distinct under deuteranopia, protanopia and tritanopia.
#   Crucially, the two-group contrast used throughout is blue #0072B2 vs
#   vermillion #D55E00 -- a high-luminance-contrast pair (not the generic
#   blue/red), so it survives grayscale printing and all three CVD types.
#
#   Diverging heatmaps use the perceptually-uniform "vik" scientific colour
#   map, which avoids the false gradients of generic blue-white-red ramps:
#     Crameri, F., Shephard, G.E., Heron, P.J. (2020) The misuse of colour in
#       science communication. Nature Communications 11:5444.
#       https://doi.org/10.1038/s41467-020-19160-7
#   vik is available exactly via scico::scale_fill_scico(palette = "vik");
#   the hard-coded ramp below samples that map so figures do not depend on
#   scico being installed.
#
#   We deliberately avoid generic software defaults (Tableau-10, ggplot2 hue,
#   jet/rainbow) and journal-brand palettes that merely re-skin the same
#   blue/red/green hues without a perceptual or accessibility basis.
# -------------------------------------------------------------------------

# Okabe-Ito Color Universal Design palette (Wong 2011, Nature Methods).
# Order: black, orange, sky blue, bluish green, yellow, blue, vermillion,
# reddish purple.
okabe_ito_palette <- function() {
  c("#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}
# Backward-compatible alias (older scripts referenced npg_palette()).
npg_palette <- function() okabe_ito_palette()

# Semantic group -> Okabe-Ito color map. Baseline/control = Okabe-Ito blue;
# perturbed/immune-high = Okabe-Ito vermillion; third group = bluish green.
group_palette <- function() {
  c(
    nucleus_pulposus = "#0072B2",  # OI blue        - control tissue
    chordoma_I       = "#D55E00",  # OI vermillion  - immune-high tumors
    chordoma_C       = "#009E73",  # OI bluish green- comparator tumors
    control          = "#0072B2",  # OI blue        - perturbation baseline
    T_DARPin         = "#D55E00",  # OI vermillion  - TBXT perturbed
    untreated        = "#0072B2",
    palbociclib      = "#D55E00",
    wild_type        = "#0072B2",
    heterozygous     = "#D55E00"
  )
}

# Diverging fill for signed-effect bar charts (delta > 0 vs delta < 0):
# Okabe-Ito vermillion (positive) vs blue (negative).
diverging_sign_fill <- function() c("TRUE" = "#D55E00", "FALSE" = "#0072B2")

# Perceptually-uniform diverging ramp sampled from Crameri's "vik"
# (Nature Communications 2020) for z-score / expression heatmaps. Use with
# ggplot2::scale_fill_gradientn(colours = vik_ramp(), ...), or take the
# endpoints from heatmap_diverging() for scale_fill_gradient2().
vik_ramp <- function() {
  c("#001261", "#033E7E", "#1E5FA3", "#6498C6", "#B3CCE0",
    "#F1E7DE", "#E0AC81", "#C87C50", "#A6492C", "#741B11", "#4D0000")
}
# Endpoints/midpoint for scale_fill_gradient2() (vik-derived).
heatmap_diverging <- function() list(low = "#1E5FA3", mid = "#F1E7DE", high = "#A6492C")

# Fill for boolean concordance/stability tiles: bluish green (concordant)
# vs vermillion (discordant) -- Okabe-Ito, CVD-safe.
concordance_fill <- function() c("TRUE" = "#009E73", "FALSE" = "#D55E00")

# -------------------------------------------------------------------------
# theme_nature(): shared publication theme for all figures.
#
# Design (Nature-family conventions):
#   - clean white background, no gridlines (data-ink maximized);
#   - thin dark axis lines and outward ticks;
#   - bold, left-aligned panel title spanning the whole plot; muted
#     subtitle/caption; compact, unobtrusive legend;
#   - sans-serif type, consistent relative sizing from a single base_size.
# Per-plot `+ theme(...)` overrides placed after this call still apply
# (e.g. blanking axes/ticks on heatmaps, moving the legend).
# -------------------------------------------------------------------------
theme_nature <- function(base_size = 10, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major   = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.background   = ggplot2::element_blank(),
      plot.background    = ggplot2::element_blank(),
      axis.line          = ggplot2::element_line(color = "grey20", linewidth = 0.3),
      axis.ticks         = ggplot2::element_line(color = "grey20", linewidth = 0.3),
      axis.ticks.length  = grid::unit(2.2, "pt"),
      axis.text          = ggplot2::element_text(color = "grey20"),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0,
                                                 size = ggplot2::rel(1.08),
                                                 margin = ggplot2::margin(b = 4)),
      plot.title.position = "plot",
      plot.subtitle      = ggplot2::element_text(color = "grey30",
                                                 size = ggplot2::rel(0.85),
                                                 margin = ggplot2::margin(b = 4)),
      plot.caption       = ggplot2::element_text(color = "grey40",
                                                 size = ggplot2::rel(0.72), hjust = 0),
      plot.caption.position = "plot",
      legend.title       = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.82)),
      legend.text        = ggplot2::element_text(size = ggplot2::rel(0.78)),
      legend.key.size    = grid::unit(10, "pt"),
      strip.text         = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      plot.margin        = ggplot2::margin(6, 8, 6, 6)
    )
}
