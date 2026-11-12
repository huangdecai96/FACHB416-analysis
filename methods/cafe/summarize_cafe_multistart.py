#!/usr/bin/env python3
"""Summarize three stochastic CAFE5 starts and select best fits per model."""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


MODELS = [
    "base_corr_l1", "error_corr_l1", "gamma3_corr_l1",
    "gamma3_error_corr_l1", "base_corr_l01", "base_relaxed_l1",
]
DATASETS = ["original", "augmented"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def find_main_results(root: Path) -> Path:
    hits = sorted(
        path for path in root.glob("*_results.txt")
        if not path.name.endswith("_family_results.txt")
        and not path.name.endswith("_clade_results.txt")
    )
    if len(hits) != 1:
        raise SystemExit(f"Expected one main *_results.txt in {root}, found {len(hits)}")
    return hits[0]


def metric(text: str, patterns: list[str]) -> str:
    for pattern in patterns:
        hit = re.search(pattern, text, flags=re.MULTILINE | re.IGNORECASE)
        if hit:
            return hit.group(1)
    return ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    root = Path(args.root)
    out = Path(args.out_dir)
    replicate_roots = {
        "replicate_01": root / "results",
        "replicate_02": root / "multistart" / "replicate_02" / "results",
        "replicate_03": root / "multistart" / "replicate_03" / "results",
    }

    fits: list[dict[str, object]] = []
    model_summaries: list[dict[str, object]] = []
    call_cache: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for replicate, replicate_root in replicate_roots.items():
        for dataset in DATASETS:
            for model in MODELS:
                model_root = replicate_root / dataset / model
                results_dir = model_root / "results"
                parsed_dir = model_root / "parsed"
                result_file = find_main_results(results_dir)
                text = result_file.read_text(encoding="utf-8")
                likelihood = metric(text, [
                    r"Likelihood \(-lnL\):\s*([0-9.eE+-]+)",
                    r"Final -lnL:\s*([0-9.eE+-]+)",
                ])
                if not likelihood:
                    raise SystemExit(f"Missing likelihood in {result_file}")
                fits.append(
                    {
                        "replicate": replicate,
                        "dataset": dataset,
                        "model": model,
                        "negative_log_likelihood": float(likelihood),
                        "lambda": metric(text, [r"^Lambda:\s*([0-9.eE+-]+)"]),
                        "alpha": metric(text, [r"^Alpha:\s*([0-9.eE+-]+)"]),
                        "result_file": result_file.relative_to(root).as_posix(),
                    }
                )
                calls = read_tsv(parsed_dir / "terminal_branch_significant_changes.tsv")
                call_cache[(replicate, dataset, model)] = calls
                for row in read_tsv(parsed_dir / "terminal_branch_summary.tsv"):
                    model_summaries.append({"replicate": replicate, "dataset": dataset, **row})

    fit_fields = [
        "replicate", "dataset", "model", "negative_log_likelihood",
        "lambda", "alpha", "result_file",
    ]
    write_tsv(out / "replicate_model_fit.tsv", fits, fit_fields)

    best: dict[tuple[str, str], dict[str, object]] = {}
    for dataset in DATASETS:
        for model in MODELS:
            candidates = [row for row in fits if row["dataset"] == dataset and row["model"] == model]
            best[(dataset, model)] = min(
                candidates,
                key=lambda row: (float(row["negative_log_likelihood"]), str(row["replicate"])),
            )
    selected_rows = [best[key] for key in sorted(best)]
    write_tsv(out / "best_fit_model_selection.tsv", selected_rows, fit_fields)

    stability_counter: dict[tuple[str, str, str, str, str], set[str]] = defaultdict(set)
    for (replicate, dataset, model), calls in call_cache.items():
        for row in calls:
            stability_counter[
                (dataset, model, row["branch_label"], row["family_id"], row["direction"])
            ].add(replicate)
    stability_rows: list[dict[str, object]] = []
    for key, replicates in sorted(stability_counter.items()):
        dataset, model, branch, family, direction = key
        stability_rows.append(
            {
                "dataset": dataset,
                "model": model,
                "branch_label": branch,
                "family_id": family,
                "direction": direction,
                "replicates_supporting": len(replicates),
                "replicates_tested": 3,
                "stable_at_least_2of3": len(replicates) >= 2,
                "supporting_replicates": ",".join(sorted(replicates)),
            }
        )
    stability_fields = [
        "dataset", "model", "branch_label", "family_id", "direction",
        "replicates_supporting", "replicates_tested", "stable_at_least_2of3",
        "supporting_replicates",
    ]
    write_tsv(out / "replicate_call_stability.tsv", stability_rows, stability_fields)

    stable_by_family: dict[tuple[str, str, str], dict[str, str]] = defaultdict(dict)
    for row in stability_rows:
        if row["stable_at_least_2of3"]:
            stable_by_family[
                (str(row["dataset"]), str(row["branch_label"]), str(row["family_id"]))
            ][str(row["model"])] = str(row["direction"])
    stable_consensus_rows: list[dict[str, object]] = []
    for (dataset, branch, family), by_model in sorted(stable_by_family.items()):
        directions = Counter(by_model.values())
        dominant, support = directions.most_common(1)[0]
        stable_consensus_rows.append(
            {
                "dataset": dataset,
                "branch_label": branch,
                "family_id": family,
                "dominant_direction": dominant,
                "stable_models_supporting_direction": support,
                "models_tested": len(MODELS),
                "replicate_stable_all_six": support == len(MODELS),
                **{
                    f"{model}_stable_direction": by_model.get(model, "not_stable_2of3")
                    for model in MODELS
                },
            }
        )
    stable_consensus_fields = [
        "dataset", "branch_label", "family_id", "dominant_direction",
        "stable_models_supporting_direction", "models_tested",
        "replicate_stable_all_six",
    ] + [f"{model}_stable_direction" for model in MODELS]
    write_tsv(
        out / "replicate_stable_six_model_family_consensus.tsv",
        stable_consensus_rows,
        stable_consensus_fields,
    )

    stable_branch_summary: list[dict[str, object]] = []
    stable_branch_labels = sorted({row["branch_label"] for row in stability_rows})
    for dataset in DATASETS:
        for branch in stable_branch_labels:
            stable_baseline = [
                row for row in stability_rows
                if row["dataset"] == dataset
                and row["model"] == "base_corr_l1"
                and row["branch_label"] == branch
                and row["stable_at_least_2of3"]
            ]
            stable_six = [
                row for row in stable_consensus_rows
                if row["dataset"] == dataset
                and row["branch_label"] == branch
                and row["replicate_stable_all_six"]
            ]
            stable_branch_summary.append(
                {
                    "dataset": dataset,
                    "branch_label": branch,
                    "stable_baseline_expanded": sum(
                        row["direction"] == "expansion" for row in stable_baseline
                    ),
                    "stable_baseline_contracted": sum(
                        row["direction"] == "contraction" for row in stable_baseline
                    ),
                    "stable_six_expanded": sum(
                        row["dominant_direction"] == "expansion" for row in stable_six
                    ),
                    "stable_six_contracted": sum(
                        row["dominant_direction"] == "contraction" for row in stable_six
                    ),
                }
            )
    stable_branch_fields = [
        "dataset", "branch_label", "stable_baseline_expanded",
        "stable_baseline_contracted", "stable_six_expanded",
        "stable_six_contracted",
    ]
    write_tsv(
        out / "replicate_stable_terminal_branch_summary.tsv",
        stable_branch_summary,
        stable_branch_fields,
    )

    best_calls: list[dict[str, object]] = []
    by_family: dict[tuple[str, str, str], dict[str, dict[str, str]]] = defaultdict(dict)
    for dataset in DATASETS:
        for model in MODELS:
            replicate = str(best[(dataset, model)]["replicate"])
            for row in call_cache[(replicate, dataset, model)]:
                selected = {
                    "dataset": dataset,
                    "selected_replicate": replicate,
                    **row,
                }
                best_calls.append(selected)
                by_family[(dataset, row["branch_label"], row["family_id"])][model] = row
    best_call_fields = ["dataset", "selected_replicate"] + list(next(iter(call_cache.values()))[0].keys())
    write_tsv(out / "best_fit_terminal_calls.tsv", best_calls, best_call_fields)

    consensus_rows: list[dict[str, object]] = []
    for (dataset, branch, family), by_model in sorted(by_family.items()):
        directions = Counter(row["direction"] for row in by_model.values())
        dominant, support = directions.most_common(1)[0]
        consensus_rows.append(
            {
                "dataset": dataset,
                "branch_label": branch,
                "family_id": family,
                "dominant_direction": dominant,
                "models_supporting_direction": support,
                "models_tested": len(MODELS),
                "robust_all_six": support == len(MODELS),
                **{
                    f"{model}_direction": by_model.get(model, {}).get("direction", "not_significant")
                    for model in MODELS
                },
            }
        )
    consensus_fields = [
        "dataset", "branch_label", "family_id", "dominant_direction",
        "models_supporting_direction", "models_tested", "robust_all_six",
    ] + [f"{model}_direction" for model in MODELS]
    write_tsv(out / "best_fit_six_model_family_consensus.tsv", consensus_rows, consensus_fields)

    branch_labels = sorted({row["branch_label"] for row in model_summaries})
    branch_summary: list[dict[str, object]] = []
    for dataset in DATASETS:
        baseline_rep = str(best[(dataset, "base_corr_l1")]["replicate"])
        baseline_calls = call_cache[(baseline_rep, dataset, "base_corr_l1")]
        for branch in branch_labels:
            baseline = [row for row in baseline_calls if row["branch_label"] == branch]
            robust = [
                row for row in consensus_rows
                if row["dataset"] == dataset and row["branch_label"] == branch and row["robust_all_six"]
            ]
            branch_summary.append(
                {
                    "dataset": dataset,
                    "branch_label": branch,
                    "baseline_selected_replicate": baseline_rep,
                    "baseline_expanded": sum(row["direction"] == "expansion" for row in baseline),
                    "baseline_contracted": sum(row["direction"] == "contraction" for row in baseline),
                    "robust_six_expanded": sum(row["dominant_direction"] == "expansion" for row in robust),
                    "robust_six_contracted": sum(row["dominant_direction"] == "contraction" for row in robust),
                }
            )
    branch_fields = [
        "dataset", "branch_label", "baseline_selected_replicate",
        "baseline_expanded", "baseline_contracted",
        "robust_six_expanded", "robust_six_contracted",
    ]
    write_tsv(out / "best_fit_terminal_branch_summary.tsv", branch_summary, branch_fields)

    model_summary_fields = ["replicate", "dataset"] + [
        "model", "branch", "branch_label", "branch_type", "significant_families",
        "expanded_families", "contracted_families", "copy_gain", "copy_loss",
    ]
    write_tsv(out / "replicate_terminal_branch_model_summary.tsv", model_summaries, model_summary_fields)

    variability_rows: list[dict[str, object]] = []
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in model_summaries:
        grouped[(str(row["dataset"]), str(row["model"]), str(row["branch_label"]))].append(row)
    for key, rows in sorted(grouped.items()):
        dataset, model, branch = key
        expanded = [int(row["expanded_families"]) for row in rows]
        contracted = [int(row["contracted_families"]) for row in rows]
        variability_rows.append(
            {
                "dataset": dataset,
                "model": model,
                "branch_label": branch,
                "expanded_min": min(expanded),
                "expanded_median": statistics.median(expanded),
                "expanded_max": max(expanded),
                "contracted_min": min(contracted),
                "contracted_median": statistics.median(contracted),
                "contracted_max": max(contracted),
            }
        )
    variability_fields = [
        "dataset", "model", "branch_label", "expanded_min", "expanded_median",
        "expanded_max", "contracted_min", "contracted_median", "contracted_max",
    ]
    write_tsv(out / "replicate_count_variability.tsv", variability_rows, variability_fields)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
