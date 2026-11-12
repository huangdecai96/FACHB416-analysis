#!/usr/bin/env python3
"""Integrate homology, translation, RNA, overlap and splice evidence for rescued models."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row[key]: row for row in rows}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-audit", type=Path, required=True)
    parser.add_argument("--translation-qc", type=Path, required=True)
    parser.add_argument("--rna-detail", type=Path, required=True)
    parser.add_argument("--overlap-detail", type=Path, required=True)
    parser.add_argument("--splice-summary", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    models = read_tsv(args.model_audit, "family_id")
    translations = read_tsv(args.translation_qc, "protein_id")
    rna = read_tsv(args.rna_detail, "protein_id")
    overlap = read_tsv(args.overlap_detail, "protein_id")
    splice = read_tsv(args.splice_summary, "protein_id")
    detail: list[dict[str, object]] = []

    for family, model in sorted(models.items()):
        protein = f"FCHB416_aug_{family}.t1"
        translation = translations.get(protein, {})
        rna_row = rna.get(protein, {})
        overlap_row = overlap.get(protein, {})
        splice_row = splice.get(protein, {})
        translation_pass = translation.get("translation_qc") == "pass"
        rna_class = rna_row.get("rna_support_class", "not_counted")
        overlap_class = overlap_row.get("annotation_overlap_class", "not_audited")
        no_substantial_overlap = overlap_class not in {
            "same_strand_CDS_overlap_ge10pct", "opposite_strand_CDS_overlap_ge10pct"
        }
        splice_class = splice_row.get("splice_support_class", "not_audited")
        splice_fraction = as_float(splice_row.get("strict_intron_support_fraction", ""), default=0.0)

        if not translation_pass or not no_substantial_overlap:
            tier = "sensitivity_only"
        elif rna_class == "strong" and splice_class == "all_introns_strictly_supported":
            tier = "high_confidence_rescue"
        elif rna_class == "strong" and splice_class == "partial_intron_support" and splice_fraction >= 0.5:
            tier = "supported_rescue_candidate"
        elif rna_class in {"strong", "moderate"}:
            tier = "provisional_rescue_candidate"
        else:
            tier = "sensitivity_only"

        detail.append(
            {
                "family_id": family,
                "protein_id": protein,
                "source_species_support": model.get("source_species_support", ""),
                "query_coverage": model.get("query_coverage", ""),
                "identity": model.get("identity", ""),
                "legacy_robust_CAFE": model.get("robust_cafe_all_models", ""),
                "translation_qc": translation.get("translation_qc", "not_counted"),
                "rna_support_class": rna_class,
                "max_CPM": rna_row.get("max_CPM", ""),
                "samples_CPM_ge_0.1": rna_row.get("samples_CPM_ge_0.1", ""),
                "annotation_overlap_class": overlap_class,
                "candidate_introns": splice_row.get("candidate_introns", ""),
                "strictly_supported_introns": splice_row.get("strictly_supported_introns", ""),
                "strict_intron_support_fraction": splice_row.get("strict_intron_support_fraction", ""),
                "splice_support_class": splice_class,
                "integrated_rescue_tier": tier,
            }
        )

    tiers = [
        "high_confidence_rescue", "supported_rescue_candidate",
        "provisional_rescue_candidate", "sensitivity_only",
    ]
    summary = [
        {"metric": "rescued_models", "value": len(detail)},
        *[
            {"metric": f"integrated_rescue_tier:{tier}", "value": sum(row["integrated_rescue_tier"] == tier for row in detail)}
            for tier in tiers
        ],
        {"metric": "legacy_robust_high_or_supported", "value": sum(
            row["legacy_robust_CAFE"] == "yes"
            and row["integrated_rescue_tier"] in {"high_confidence_rescue", "supported_rescue_candidate"}
            for row in detail
        )},
    ]
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "candidate_rescue_integrated_evidence.tsv", detail, list(detail[0]))
    write_tsv(args.outdir / "candidate_rescue_integrated_summary.tsv", summary, ["metric", "value"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
