#!/usr/bin/env python3
"""Parse significant CAFE5 changes for every branch, including all terminals."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def find_one(results: Path, suffix: str) -> Path:
    hits = sorted(results.glob(f"*{suffix}"))
    if len(hits) != 1:
        raise SystemExit(f"Expected one *{suffix} in {results}, found {len(hits)}")
    return hits[0]


def optional_float(value: str | None) -> float:
    if value in (None, "", "N/A", "nan"):
        return math.nan
    return float(value)


def branch_type(branch: str) -> tuple[str, str]:
    match = re.fullmatch(r"(.*)<(\d+)>", branch)
    if not match:
        return "unknown", branch
    label = match.group(1)
    return ("terminal", label) if label else ("internal", branch)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-name", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--pvalue", type=float, default=0.05)
    args = parser.parse_args()

    results = Path(args.results_dir)
    out = Path(args.out_dir)
    family_path = find_one(results, "_family_results.txt")
    prefix = family_path.name.removesuffix("_family_results.txt")
    family_rows = read_tsv(family_path)
    family_p = {row["#FamilyID"]: float(row["pvalue"]) for row in family_rows}
    globally_sig = {
        row["#FamilyID"]
        for row in family_rows
        if row.get("Significant at 0.05", "").lower() == "y"
        and float(row["pvalue"]) <= args.pvalue
    }

    changes = read_tsv(results / f"{prefix}_change.tab")
    probs = read_tsv(results / f"{prefix}_branch_probabilities.tab")
    if not changes:
        raise SystemExit("No CAFE change rows")
    branches = [column for column in changes[0] if column != "FamilyID"]
    prob_by_family = {row["FamilyID"]: row for row in probs}

    calls: list[dict[str, object]] = []
    for row in changes:
        family = row["FamilyID"]
        if family not in globally_sig:
            continue
        for branch in branches:
            change = int(row[branch])
            probability = optional_float(prob_by_family.get(family, {}).get(branch))
            if change == 0 or math.isnan(probability) or probability > args.pvalue:
                continue
            kind, label = branch_type(branch)
            calls.append(
                {
                    "model": args.model_name,
                    "cafe_prefix": prefix,
                    "family_id": family,
                    "branch": branch,
                    "branch_label": label,
                    "branch_type": kind,
                    "change": change,
                    "direction": "expansion" if change > 0 else "contraction",
                    "branch_pvalue": probability,
                    "family_pvalue": family_p[family],
                }
            )

    fields = [
        "model", "cafe_prefix", "family_id", "branch", "branch_label", "branch_type",
        "change", "direction", "branch_pvalue", "family_pvalue",
    ]
    write_tsv(out / "all_branch_significant_changes.tsv", calls, fields)
    terminal = [row for row in calls if row["branch_type"] == "terminal"]
    write_tsv(out / "terminal_branch_significant_changes.tsv", terminal, fields)

    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in calls:
        grouped[(str(row["branch"]), str(row["branch_label"]), str(row["branch_type"]))].append(row)
    summary: list[dict[str, object]] = []
    for branch in branches:
        kind, label = branch_type(branch)
        rows = grouped.get((branch, label, kind), [])
        summary.append(
            {
                "model": args.model_name,
                "branch": branch,
                "branch_label": label,
                "branch_type": kind,
                "significant_families": len(rows),
                "expanded_families": sum(row["direction"] == "expansion" for row in rows),
                "contracted_families": sum(row["direction"] == "contraction" for row in rows),
                "copy_gain": sum(max(0, int(row["change"])) for row in rows),
                "copy_loss": sum(min(0, int(row["change"])) for row in rows),
            }
        )
    summary_fields = [
        "model", "branch", "branch_label", "branch_type", "significant_families",
        "expanded_families", "contracted_families", "copy_gain", "copy_loss",
    ]
    write_tsv(out / "all_branch_summary.tsv", summary, summary_fields)
    write_tsv(out / "terminal_branch_summary.tsv",
              [row for row in summary if row["branch_type"] == "terminal"], summary_fields)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
