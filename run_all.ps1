$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

$steps = @(
    "scripts/01_gse230168_pathway_scores.R",
    "scripts/02_gse230168_microenvironment_biomarker_context.R",
    "scripts/03_gse276715_tbxt_signature.R",
    "scripts/04_gse216418_palbociclib_signature.R",
    "scripts/05_gse312699_naive_validation.R",
    "scripts/06_gse239531_confirmatory_context.R",
    "scripts/07_yin2024_proteomic_input_qc.R",
    "scripts/08_yin2024_unsupervised_crossmodality.R",
    "scripts/09_yin2024_supervised_labeltransfer.R",
    "scripts/10_gse277707_developmental_context.R",
    "scripts/11_signature_robustness_and_biomarkers.R",
    "scripts/12_effect_size_stability_analysis.R",
    "scripts/13_external_evidence_and_druggability.R",
    "scripts/14_analysis_figures_and_tables.R",
    "scripts/15_validation_qc_outputs.R"
)

foreach ($script in $steps) {
    Write-Host "Running $script"
    & Rscript $script
    if ($LASTEXITCODE -ne 0) {
        throw "Pipeline step failed with exit code ${LASTEXITCODE}: $script"
    }
}

Write-Host "Analysis pipeline completed."
