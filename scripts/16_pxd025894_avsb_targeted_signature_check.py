#!/usr/bin/env python3
"""
Targeted PXD025894/Shen 2021 AvsB support check.

D16 scope: use the published 258-protein AvsB differential supplement only as
supporting targeted proteomic evidence. This script reports row-level Shen
A/B direction and p-values for AvsB proteins that are present in the selected
Hallmark/custom signature universe used by scripts/00_chordoma_utils.R.
It does not score absent proteins and does not make a proteome-wide claim.
"""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
AVSB_CSV = ROOT / "data" / "processed" / "PXD025894_Shen2021" / "Shen2021_AvsB_differential_proteins_wide.csv"
CUSTOM_CSV = ROOT / "results" / "custom_signature_genes.csv"
HALLMARK_SYMBOLS_GMT = ROOT / "data" / "raw" / "h.all.v2025.1.Hs.symbols.gmt"

OUT_HITS = ROOT / "results" / "PXD025894_AvsB_targeted_signature_hits.csv"
OUT_SUMMARY = ROOT / "results" / "PXD025894_AvsB_targeted_signature_summary.csv"
OUT_REPORT = ROOT / "results" / "PXD025894_AvsB_targeted_signature_report.md"

SELECTED_HALLMARKS = [
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
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
]

CUSTOM_SIGNATURES = [
    "CUSTOM_ANTIGEN_PRESENTATION",
    "CUSTOM_CHORDOMA_MARKERS",
    "CUSTOM_CDK46_AXIS",
    "CUSTOM_PLK1_AXIS",
    "CUSTOM_DRUGGABLE_RTKS",
]

SIGNATURE_LABELS = {
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "Interferon alpha",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "Interferon gamma",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": "IL6-JAK-STAT3",
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
    "HALLMARK_COMPLEMENT": "Complement",
    "HALLMARK_E2F_TARGETS": "E2F targets",
    "HALLMARK_G2M_CHECKPOINT": "G2M checkpoint",
    "HALLMARK_MTORC1_SIGNALING": "mTORC1 signaling",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
    "HALLMARK_ANGIOGENESIS": "Angiogenesis",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "WNT/beta-catenin",
    "CUSTOM_ANTIGEN_PRESENTATION": "Antigen presentation",
    "CUSTOM_CHORDOMA_MARKERS": "Chordoma markers",
    "CUSTOM_CDK46_AXIS": "CDK4/6 axis",
    "CUSTOM_PLK1_AXIS": "PLK1 axis",
    "CUSTOM_DRUGGABLE_RTKS": "Druggable RTKs",
}

SIGNATURE_AXIS = {
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "immune/interferon",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "immune/interferon",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": "immune/interferon",
    "HALLMARK_INFLAMMATORY_RESPONSE": "immune/interferon",
    "HALLMARK_COMPLEMENT": "immune/interferon",
    "HALLMARK_E2F_TARGETS": "cell-cycle",
    "HALLMARK_G2M_CHECKPOINT": "cell-cycle",
    "HALLMARK_MTORC1_SIGNALING": "metabolic/stress context",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT/stromal context",
    "HALLMARK_ANGIOGENESIS": "RTK/angiogenic context",
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING": "developmental context",
    "CUSTOM_ANTIGEN_PRESENTATION": "antigen presentation / immune state",
    "CUSTOM_CHORDOMA_MARKERS": "chordoma identity",
    "CUSTOM_CDK46_AXIS": "cell-cycle / CDK4-6",
    "CUSTOM_PLK1_AXIS": "cell-cycle / PLK1",
    "CUSTOM_DRUGGABLE_RTKS": "RTK/angiogenic context",
}


def normalize_symbol(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().upper()


def parse_gene_symbol(row: pd.Series) -> str:
    symbol = normalize_symbol(row.get("Gene name"))
    if symbol:
        return symbol
    desc = str(row.get("Protein description", ""))
    match = re.search(r"\bGN=([^\s]+)", desc)
    return normalize_symbol(match.group(1)) if match else ""


def read_hallmark_symbols() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with HALLMARK_SYMBOLS_GMT.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            if name in SELECTED_HALLMARKS:
                out[name] = {normalize_symbol(g) for g in parts[2:] if normalize_symbol(g)}
    missing = sorted(set(SELECTED_HALLMARKS) - set(out))
    if missing:
        raise RuntimeError(f"Selected Hallmark sets missing from GMT: {missing}")
    return out


def read_custom_signatures() -> dict[str, set[str]]:
    df = pd.read_csv(CUSTOM_CSV, dtype=str)
    df = df[df["signature_name"].isin(CUSTOM_SIGNATURES)].copy()
    out: dict[str, set[str]] = {}
    for name, sub in df.groupby("signature_name", sort=False):
        out[name] = {normalize_symbol(g) for g in sub["gene_symbol"] if normalize_symbol(g)}
    missing = sorted(set(CUSTOM_SIGNATURES) - set(out))
    if missing:
        raise RuntimeError(f"Selected custom signatures missing from table: {missing}")
    return out


def load_avsb() -> pd.DataFrame:
    df = pd.read_csv(AVSB_CSV, dtype=str)
    required = ["Protein accession", "Protein description", "A/B Ratio", "Regulated Type", "A/B P value"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"AvsB table missing required columns: {missing}")
    df["gene_symbol"] = df.apply(parse_gene_symbol, axis=1)
    df["ab_ratio"] = pd.to_numeric(df["A/B Ratio"], errors="coerce")
    df["log2_ab_ratio"] = df["ab_ratio"].map(lambda x: math.log2(x) if pd.notna(x) and x > 0 else float("nan"))
    df["p_value"] = pd.to_numeric(df["A/B P value"], errors="coerce")
    df["direction"] = df["Regulated Type"].astype(str).str.strip()
    return df


def format_gene_list(values: list[str], max_items: int = 20) -> str:
    values = [v for v in values if v]
    if len(values) <= max_items:
        return "; ".join(values)
    return "; ".join(values[:max_items]) + f"; ... (+{len(values) - max_items} more)"


def main() -> None:
    OUT_HITS.parent.mkdir(parents=True, exist_ok=True)

    signature_sets = {}
    signature_sets.update(read_hallmark_symbols())
    signature_sets.update(read_custom_signatures())
    signature_sets = {k: signature_sets[k] for k in SELECTED_HALLMARKS + CUSTOM_SIGNATURES}

    avsb = load_avsb()
    avsb_gene_universe = set(avsb["gene_symbol"]) - {""}

    records: list[dict[str, object]] = []
    for signature_name, genes in signature_sets.items():
        present = sorted(genes & avsb_gene_universe)
        if not present:
            continue
        sub = avsb[avsb["gene_symbol"].isin(present)].copy()
        sub = sub.sort_values(["p_value", "gene_symbol", "Protein accession"], na_position="last")
        for _, row in sub.iterrows():
            records.append(
                {
                    "signature_name": signature_name,
                    "signature_label": SIGNATURE_LABELS.get(signature_name, signature_name),
                    "axis_or_context": SIGNATURE_AXIS.get(signature_name, ""),
                    "signature_family": "selected Hallmark" if signature_name.startswith("HALLMARK_") else "custom pathway/signature",
                    "signature_genes_total": len(genes),
                    "signature_genes_present_in_avsb": len(present),
                    "gene_symbol": row["gene_symbol"],
                    "protein_accession": row["Protein accession"],
                    "protein_description": row["Protein description"],
                    "ab_ratio": row["ab_ratio"],
                    "log2_ab_ratio": row["log2_ab_ratio"],
                    "direction": row["direction"],
                    "p_value": row["p_value"],
                    "subcellular_localization": row.get("Subcellular localization", ""),
                }
            )

    hits = pd.DataFrame.from_records(records)
    if hits.empty:
        raise RuntimeError("No selected signature genes were found in AvsB table.")

    hits.to_csv(OUT_HITS, index=False, quoting=csv.QUOTE_MINIMAL)

    summary_records = []
    for signature_name, genes in signature_sets.items():
        sub = hits[hits["signature_name"] == signature_name].copy()
        present_genes = sorted(set(sub["gene_symbol"])) if not sub.empty else []
        up_genes = sorted(set(sub.loc[sub["direction"].str.upper() == "UP", "gene_symbol"])) if not sub.empty else []
        down_genes = sorted(set(sub.loc[sub["direction"].str.upper() == "DOWN", "gene_symbol"])) if not sub.empty else []
        sig_p = sub[sub["p_value"] < 0.05].copy() if not sub.empty else sub
        summary_records.append(
            {
                "signature_name": signature_name,
                "signature_label": SIGNATURE_LABELS.get(signature_name, signature_name),
                "axis_or_context": SIGNATURE_AXIS.get(signature_name, ""),
                "signature_family": "selected Hallmark" if signature_name.startswith("HALLMARK_") else "custom pathway/signature",
                "signature_genes_total": len(genes),
                "n_genes_present_in_avsb": len(present_genes),
                "coverage_fraction": len(present_genes) / len(genes) if genes else 0,
                "n_up": len(up_genes),
                "n_down": len(down_genes),
                "n_p_lt_0_05": len(set(sig_p["gene_symbol"])) if not sig_p.empty else 0,
                "min_p_value": sub["p_value"].min() if not sub.empty else float("nan"),
                "present_genes": ";".join(present_genes),
                "up_genes": ";".join(up_genes),
                "down_genes": ";".join(down_genes),
            }
        )
    summary = pd.DataFrame.from_records(summary_records)
    summary.to_csv(OUT_SUMMARY, index=False, quoting=csv.QUOTE_MINIMAL)

    unique_hit_genes = sorted(set(hits["gene_symbol"]))
    asns_rows = hits[hits["gene_symbol"] == "ASNS"].copy()
    asns_avsb = avsb[avsb["gene_symbol"] == "ASNS"].copy()
    top_rows = hits.sort_values(["p_value", "signature_name", "gene_symbol"], na_position="last").head(20)
    custom_rows = hits[hits["signature_family"] == "custom pathway/signature"].copy()

    report = []
    report.append("# PXD025894/Shen 2021 AvsB targeted signature check")
    report.append("")
    report.append("## Scope")
    report.append("")
    report.append("D16 scope: PXD025894/Shen 2021 is supporting targeted proteomic evidence only. This report uses the published 258-protein AvsB differential supplement and reports Shen row-level A/B direction and p-values only for proteins that are present in the selected signature universe. It does not impute absent proteins, compute proteome-wide scores, or support a primary validation claim.")
    report.append("")
    report.append("## Inputs")
    report.append("")
    report.append(f"- AvsB table: `{AVSB_CSV.relative_to(ROOT)}`")
    report.append(f"- Selected Hallmark signatures: `{HALLMARK_SYMBOLS_GMT.relative_to(ROOT)}`")
    report.append(f"- Custom signatures: `{CUSTOM_CSV.relative_to(ROOT)}`")
    report.append("")
    report.append("## Basic checks")
    report.append("")
    report.append(f"- AvsB rows: {len(avsb)}")
    report.append(f"- AvsB rows with parsed gene symbol: {int((avsb['gene_symbol'] != '').sum())}")
    report.append(f"- Selected signatures checked: {len(signature_sets)} (11 Hallmark + 5 custom)")
    report.append(f"- Unique selected signature genes: {len(set().union(*signature_sets.values()))}")
    report.append(f"- Unique AvsB genes intersecting selected signatures: {len(unique_hit_genes)}")
    report.append(f"- Long-form signature membership hits: {len(hits)}")
    report.append("")
    report.append("## ASNS positive-control check")
    report.append("")
    if not asns_avsb.empty:
        row = asns_avsb.iloc[0]
        memberships = sorted(set(asns_rows["signature_name"])) if not asns_rows.empty else []
        report.append(
            f"ASNS is present in AvsB: A/B ratio {row['ab_ratio']:.3f}, direction {row['direction']}, p={row['p_value']:.6g}. Selected-signature memberships found here: {', '.join(memberships) if memberships else 'none in selected signatures'}."
        )
    else:
        report.append("ASNS was not found in the parsed AvsB table; this should be investigated before using the supplement as a targeted check.")
    report.append("")
    report.append("## Per-signature coverage")
    report.append("")
    report.append("| Signature | Present / total | Up | Down | p<0.05 | Genes present |")
    report.append("|---|---:|---:|---:|---:|---|")
    for _, row in summary.iterrows():
        report.append(
            f"| {row['signature_name']} | {int(row['n_genes_present_in_avsb'])}/{int(row['signature_genes_total'])} | {int(row['n_up'])} | {int(row['n_down'])} | {int(row['n_p_lt_0_05'])} | {format_gene_list(str(row['present_genes']).split(';'))} |"
        )
    report.append("")
    report.append("## Most significant row-level hits")
    report.append("")
    report.append("These are not adjusted enrichment tests; they are Shen supplement row-level p-values for signature-member proteins present in AvsB.")
    report.append("")
    report.append("| Signature | Gene | A/B ratio | Direction | p-value |")
    report.append("|---|---:|---:|---|---:|")
    for _, row in top_rows.iterrows():
        report.append(f"| {row['signature_name']} | {row['gene_symbol']} | {row['ab_ratio']:.3f} | {row['direction']} | {row['p_value']:.6g} |")
    report.append("")
    report.append("## Custom signature hits")
    report.append("")
    if custom_rows.empty:
        report.append("No custom-signature genes were present in AvsB.")
    else:
        report.append("| Signature | Gene | A/B ratio | Direction | p-value |")
        report.append("|---|---:|---:|---|---:|")
        for _, row in custom_rows.sort_values(["signature_name", "p_value", "gene_symbol"], na_position="last").iterrows():
            report.append(f"| {row['signature_name']} | {row['gene_symbol']} | {row['ab_ratio']:.3f} | {row['direction']} | {row['p_value']:.6g} |")
    report.append("")
    report.append("## Outputs")
    report.append("")
    report.append(f"- Long-form hits: `{OUT_HITS.relative_to(ROOT)}`")
    report.append(f"- Per-signature summary: `{OUT_SUMMARY.relative_to(ROOT)}`")
    report.append(f"- Report: `{OUT_REPORT.relative_to(ROOT)}`")
    report.append("")

    OUT_REPORT.write_text("\n".join(report), encoding="utf-8")

    print(f"Wrote {OUT_HITS.relative_to(ROOT)}")
    print(f"Wrote {OUT_SUMMARY.relative_to(ROOT)}")
    print(f"Wrote {OUT_REPORT.relative_to(ROOT)}")
    print(f"AvsB rows: {len(avsb)}")
    print(f"Unique AvsB genes intersecting selected signatures: {len(unique_hit_genes)}")
    if not asns_avsb.empty:
        row = asns_avsb.iloc[0]
        print(f"ASNS: A/B={row['ab_ratio']:.3f}; direction={row['direction']}; p={row['p_value']:.6g}")


if __name__ == "__main__":
    main()