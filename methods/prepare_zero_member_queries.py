#!/usr/bin/env python3
"""Select close-species representatives for contracted families lacking FCHB416 members."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


PREFERRED = ["DOE0152z", "UTEX393", "NREL46BD3", "Mneglectum"]


def fasta_dict(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    ident = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if ident is not None:
                    result[ident] = "".join(chunks).rstrip("*")
                ident = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if ident is not None:
        result[ident] = "".join(chunks).rstrip("*")
    return result


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contracted", required=True)
    ap.add_argument("--orthogroups", required=True)
    ap.add_argument("--proteome-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    contracted = {r["family_id"]: r for r in read_tsv(Path(args.contracted))}
    og_rows = {r["Orthogroup"]: r for r in read_tsv(Path(args.orthogroups)) if r["Orthogroup"] in contracted}
    species = [x for x in next(iter(og_rows.values())) if x != "Orthogroup"]
    proteomes = {s: fasta_dict(Path(args.proteome_dir) / f"{s}.fa") for s in species}

    family_rows = []
    query_rows = []
    query_fasta: list[str] = []
    for family, cmeta in sorted(contracted.items()):
        row = og_rows[family]
        if row.get("FCHB416", "").strip():
            continue
        selected: list[tuple[str, str, str]] = []
        for sp in PREFERRED:
            members = [x.strip() for x in row.get(sp, "").split(",") if x.strip()]
            available = [(m, proteomes[sp].get(m, "")) for m in members]
            available = [(m, seq) for m, seq in available if seq]
            if available:
                member, seq = max(available, key=lambda x: len(x[1]))
                selected.append((sp, member, seq))
        if not selected:
            fallback = []
            for sp in species:
                if sp == "FCHB416":
                    continue
                for member in [x.strip() for x in row.get(sp, "").split(",") if x.strip()]:
                    seq = proteomes[sp].get(member, "")
                    if seq:
                        fallback.append((sp, member, seq))
            if fallback:
                selected = [max(fallback, key=lambda x: len(x[2]))]
        all_member_count = sum(len([x for x in row.get(sp, "").split(",") if x.strip()]) for sp in species)
        family_rows.append(
            {
                "family_id": family,
                "cafe_change": cmeta["change"],
                "total_non_fchb_members": all_member_count,
                "selected_queries": len(selected),
                "selected_species": ",".join(x[0] for x in selected),
            }
        )
        for index, (sp, member, seq) in enumerate(selected, 1):
            qid = f"{family}__{sp}__q{index}"
            query_rows.append(
                {
                    "query_id": qid,
                    "family_id": family,
                    "source_species": sp,
                    "source_protein_id": member,
                    "protein_length": len(seq),
                }
            )
            query_fasta.extend([f">{qid}", seq])

    for name, rows in [("zero_member_families.tsv", family_rows), ("zero_member_query_map.tsv", query_rows)]:
        with (out / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    (out / "zero_member_queries.faa").write_text("\n".join(query_fasta) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
