# Reproducibility guide

## Environment

The pipeline was developed with R 4.3.3 on Windows/PowerShell. Core R dependencies are:

- ggplot2
- dplyr
- readr
- tidyr
- tibble
- cowplot
- ggrepel
- scales
- readxl
- limma from Bioconductor

The optional Python script requires pandas; install from requirements.txt.

No controlled-access dataset is required.

## Public input paths

Place these files under data/raw/:

- GSE230168_count_table.csv
- GSE230168_spl_clust.csv
- GSE276715/GSE276715_RAW.tar
- GSE216418/GSE216418_RAW.tar
- GSE312699/GSE312699_series_matrix.txt.gz
- GSE312699/GPL15207_family.soft.gz
- GSE239531_raw_counts.tsv.gz (the script attempts GEO download if absent)
- GSE277707/GSE277707_norm_matrix_48and72hrs.csv.gz
- Yin2024_protein_abundance.csv
- Yin2024_mmc2_clinical_cluster.xlsx
- Yin2024_mmc4_cluster_markers.xlsx
- h.all.v2025.1.Hs.entrez.gmt
- NCBI_gene2ensembl.gz

Depending on local identifier mapping, hgnc_complete_set.tsv may be used as a fallback for GSE239531.

For the optional PXD025894 targeted check, also provide:

- data/raw/h.all.v2025.1.Hs.symbols.gmt
- data/processed/PXD025894_Shen2021/Shen2021_AvsB_differential_proteins_wide.csv

The AvsB file is the published 258-protein differential supplement, not a complete proteome-wide abundance matrix.

## Ordered analysis

run_all.ps1 and run_all.sh run scripts 01 through 15 and stop at the first error. Script 00 contains shared functions.

1. 01_gse230168_pathway_scores.R
2. 02_gse230168_microenvironment_biomarker_context.R
3. 03_gse276715_tbxt_signature.R
4. 04_gse216418_palbociclib_signature.R
5. 05_gse312699_naive_validation.R
6. 06_gse239531_confirmatory_context.R
7. 07_yin2024_proteomic_input_qc.R
8. 08_yin2024_unsupervised_crossmodality.R
9. 09_yin2024_supervised_labeltransfer.R
10. 10_gse277707_developmental_context.R
11. 11_signature_robustness_and_biomarkers.R
12. 12_effect_size_stability_analysis.R
13. 13_external_evidence_and_druggability.R
14. 14_analysis_figures_and_tables.R
15. 15_validation_qc_outputs.R

The optional 16_pxd025894_avsb_targeted_signature_check.py is not in the default runner because it requires a separately prepared published supplement.

## Statistical conventions

- Gene identifiers are harmonized to Entrez IDs.
- Count matrices use low-expression filtering followed by log2(CPM + 1).
- TPM or normalized matrices use log2(value + 1) where applicable.
- GSE312699 uses paired primary-versus-recurrent tests by patient.
- Pathway scores are mean sample-wise z-scores.
- Sensitivity analyses use within-sample percentile ranks.
- Group tests use Wilcoxon tests with Benjamini-Hochberg adjustment.
- Gene-level GSE276715, GSE216418, and GSE312699 signatures use limma.
- Effect-size stability uses Cliff's delta, bootstrap confidence intervals, and leave-one-out direction checks.
- Cross-dataset interpretation occurs at signature/pathway level without batch-merging matrices.
- Bulk microenvironment scores are marker-set summaries, not cell-type deconvolution.

## Evidence roles

GSE230168 is the patient-tumor anchor. GSE312699 is the held-out paired primary/recurrent RNA validation cohort. GSE239531 and Yin 2024 are confirmatory context layers. GSE276715 and GSE216418 provide functional perturbation evidence. GSE277707 provides developmental context. PXD025894/Shen 2021 is an optional targeted proteomic check limited to proteins present in the published AvsB differential table.

## Generated content

Scripts write intermediate matrices to data/processed/, tables and reports to results/, and graphics to figures/. These directories are ignored because their contents are regenerated locally.
