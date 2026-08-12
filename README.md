# Chordoma candidate states

Reproducible public-data analysis code for the study **An evidence-layer framework for public-omics integration in rare cancers:immune/interferon and TBXT-linked candidate states in chordoma**.

The repository integrates public transcriptomic and proteomic datasets at the signature and pathway level. It contains analysis code only: no manuscript generator, cover-letter generator, submission packaging, editorial workflow, or internal task scripts are included.

## Repository structure

- scripts/: ordered R analysis pipeline plus one optional Python targeted proteomic check.
- data/: local public inputs and generated intermediate matrices; only its README is versioned.
- results/: regenerated numerical outputs; only its README is versioned.
- figures/: regenerated figures; only its README is versioned.
- README_reproducibility.md: input paths, dependencies, execution order, and interpretation boundaries.

## Quick start

1. Install R 4.3 or later and the packages listed in README_reproducibility.md.
2. Place the public input files under data/raw/ using the documented names.
3. Run the ordered pipeline from the repository root:

    pwsh ./run_all.ps1

On macOS or Linux:

    bash ./run_all.sh

The optional PXD025894/Shen 2021 targeted check is run separately after its processed AvsB input and symbol-level Hallmark GMT are available:

    python scripts/16_pxd025894_avsb_targeted_signature_check.py

## Analysis scope

The pipeline covers:

- patient-tumor discovery in GSE230168;
- TBXT/T-DARPin perturbation in GSE276715;
- palbociclib PDX perturbation in GSE216418;
- held-out paired primary/recurrent validation in GSE312699;
- confirmatory RNA context in GSE239531;
- Yin 2024 proteomic input QC, unsupervised analysis, and marker-based label transfer;
- developmental context in GSE277707;
- rank-based sensitivity, effect-size stability, biomarker, external-evidence, and validation-QC outputs.

## Interpretation boundary

Cross-study evidence is integrated at the signature/pathway and direction-of-effect level; expression matrices are not batch-merged. Results are hypothesis-generating and do not establish clinical treatment response or a validated clinical classifier.

All required study inputs are public and de-identified. Large source files, generated matrices, result tables, and figures are intentionally not versioned.

## License

MIT License. See LICENSE.
