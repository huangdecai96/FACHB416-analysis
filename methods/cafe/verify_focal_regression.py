#!/usr/bin/env python3
"""Audit old exact counts while requiring contraction-dominated reproduction."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


EXPECTED = {
    "original": {
        "baseline_expanded": 89,
        "baseline_contracted": 248,
        "robust_six_expanded": 42,
        "robust_six_contracted": 131,
    },
    "augmented": {
        "baseline_expanded": 90,
        "baseline_contracted": 233,
        "robust_six_expanded": 28,
        "robust_six_contracted": 97,
    },
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--results-root")
    parser.add_argument("--strict-exact", action="store_true")
    args = parser.parse_args()
    root = Path(args.root)
    results_root = Path(args.results_root) if args.results_root else root / "results"
    audit_rows: list[dict[str, object]] = []
    failures: list[str] = []
    exact_failures: list[str] = []

    for dataset, expected in EXPECTED.items():
        path = results_root / dataset / "summary" / "terminal_branch_consensus_summary.tsv"
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        focal = next((row for row in rows if row["branch_label"] == "FCHB416"), None)
        if focal is None:
            failures.append(f"{dataset}: FCHB416 terminal branch missing")
            continue
        for metric, expected_value in expected.items():
            observed = int(focal[metric])
            passed = observed == expected_value
            audit_rows.append(
                {
                    "dataset": dataset,
                    "branch": "FCHB416",
                    "metric": metric,
                    "expected": expected_value,
                    "observed": observed,
                    "status": "PASS" if passed else "FAIL",
                }
            )
            if not passed:
                exact_failures.append(
                    f"{dataset}/{metric}: expected {expected_value}, observed {observed}"
                )
        baseline_directional = (
            int(focal["baseline_contracted"]) > int(focal["baseline_expanded"])
        )
        robust_directional = (
            int(focal["robust_six_contracted"]) > int(focal["robust_six_expanded"])
            and int(focal["robust_six_total"]) > 0
        )
        if not baseline_directional or not robust_directional:
            failures.append(
                f"{dataset}: contraction-dominated direction was not reproduced"
            )

    audit_path = results_root / "focal_regression_audit.tsv"
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["dataset", "branch", "metric", "expected", "observed", "status"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)

    if args.strict_exact and exact_failures:
        failures.extend(exact_failures)
    if failures:
        raise SystemExit("Formal-result regression failed: " + "; ".join(failures))
    print("FCHB416 contraction-dominated directional regression: PASS")
    if exact_failures:
        print("Exact old-count reproduction: NO")
        print("; ".join(exact_failures))
    else:
        print("Exact old-count reproduction: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
