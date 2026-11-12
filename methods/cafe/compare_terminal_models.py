#!/usr/bin/env python3
"""Build all-model consensus summaries for every terminal CAFE5 branch."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


MODELS = [
    "base_corr_l1", "error_corr_l1", "gamma3_corr_l1",
    "gamma3_error_corr_l1", "base_corr_l01", "base_relaxed_l1",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--baseline", default="base_corr_l1")
    args = parser.parse_args()

    root = Path(args.dataset_dir)
    out = Path(args.out_dir)
    calls: dict[tuple[str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    model_summaries: list[dict[str, str]] = []
    for model in MODELS:
        call_path = root / model / "parsed" / "terminal_branch_significant_changes.tsv"
        summary_path = root / model / "parsed" / "terminal_branch_summary.tsv"
        if not call_path.exists() or not summary_path.exists():
            raise SystemExit(f"Missing parsed output for model {model}")
        for row in read_tsv(call_path):
            calls[(row["branch_label"], row["family_id"])][model] = row
        model_summaries.extend(read_tsv(summary_path))

    consensus: list[dict[str, object]] = []
    for (branch, family), by_model in sorted(calls.items()):
        directions = Counter(row["direction"] for row in by_model.values())
        dominant, support = directions.most_common(1)[0]
        baseline_direction = by_model.get(args.baseline, {}).get("direction", "not_significant")
        consensus.append(
            {
                "branch_label": branch,
                "family_id": family,
                "dominant_direction": dominant,
                "models_supporting_direction": support,
                "models_with_significant_call": len(by_model),
                "models_tested": len(MODELS),
                "support_fraction": support / len(MODELS),
                "robust_all_six": support == len(MODELS),
                "baseline_direction": baseline_direction,
                **{
                    f"{model}_direction": by_model.get(model, {}).get("direction", "not_significant")
                    for model in MODELS
                },
                **{
                    f"{model}_change": by_model.get(model, {}).get("change", "0")
                    for model in MODELS
                },
            }
        )

    consensus_fields = [
        "branch_label", "family_id", "dominant_direction", "models_supporting_direction",
        "models_with_significant_call", "models_tested", "support_fraction",
        "robust_all_six", "baseline_direction",
    ] + [f"{model}_direction" for model in MODELS] + [f"{model}_change" for model in MODELS]
    write_tsv(out / "terminal_family_consensus.tsv", consensus, consensus_fields)

    robust_summary: list[dict[str, object]] = []
    branches = sorted({row["branch_label"] for row in model_summaries})
    for branch in branches:
        rows = [row for row in consensus if row["branch_label"] == branch]
        robust = [row for row in rows if row["robust_all_six"]]
        baseline = [row for row in rows if row["baseline_direction"] != "not_significant"]
        robust_summary.append(
            {
                "branch_label": branch,
                "baseline_expanded": sum(row["baseline_direction"] == "expansion" for row in baseline),
                "baseline_contracted": sum(row["baseline_direction"] == "contraction" for row in baseline),
                "robust_six_expanded": sum(row["dominant_direction"] == "expansion" for row in robust),
                "robust_six_contracted": sum(row["dominant_direction"] == "contraction" for row in robust),
                "robust_six_total": len(robust),
            }
        )
    write_tsv(
        out / "terminal_branch_consensus_summary.tsv",
        robust_summary,
        ["branch_label", "baseline_expanded", "baseline_contracted",
         "robust_six_expanded", "robust_six_contracted", "robust_six_total"],
    )
    write_tsv(
        out / "terminal_branch_model_summary.tsv",
        model_summaries,
        ["model", "branch", "branch_label", "branch_type", "significant_families",
         "expanded_families", "contracted_families", "copy_gain", "copy_loss"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
