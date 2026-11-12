#!/usr/bin/env python3
"""Concatenate trimmed single-copy protein alignments with strict taxon checks."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = ""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current:
                records[current].append(line.strip())
    return {name: "".join(parts) for name, parts in records.items()}


def write_fasta(path: Path, records: dict[str, str], width: int = 80) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start:start + width] + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--partition-tsv", type=Path, required=True)
    parser.add_argument("--summary-tsv", type=Path, required=True)
    args = parser.parse_args()

    files = sorted(args.alignment_dir.glob("*.trim.fa"))
    if not files:
        raise SystemExit(f"No *.trim.fa files in {args.alignment_dir}")

    species_order: list[str] = []
    concatenated: dict[str, list[str]] = {}
    partitions: list[dict[str, object]] = []
    cursor = 1
    for index, path in enumerate(files):
        raw = parse_fasta(path)
        by_species: dict[str, str] = {}
        for header, sequence in raw.items():
            species = header.split("__", 1)[0]
            if species in by_species:
                raise SystemExit(f"Duplicate species {species} in {path}")
            by_species[species] = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY\-X]", "X", sequence.upper())
        lengths = {len(sequence) for sequence in by_species.values()}
        if len(lengths) != 1:
            raise SystemExit(f"Unequal aligned lengths in {path}: {sorted(lengths)}")
        length = lengths.pop()
        if length == 0:
            raise SystemExit(f"Empty trimmed alignment: {path}")
        if index == 0:
            species_order = sorted(by_species)
            concatenated = {species: [] for species in species_order}
        if set(by_species) != set(species_order):
            missing = sorted(set(species_order) - set(by_species))
            extra = sorted(set(by_species) - set(species_order))
            raise SystemExit(f"Taxon mismatch in {path}: missing={missing}, extra={extra}")
        for species in species_order:
            concatenated[species].append(by_species[species])
        family = path.name.replace(".trim.fa", "")
        partitions.append({"orthogroup": family, "start_aa": cursor, "end_aa": cursor + length - 1, "length_aa": length})
        cursor += length

    joined = {species: "".join(parts) for species, parts in concatenated.items()}
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(args.output_fasta, joined)
    with args.partition_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(partitions[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(partitions)
    summary = [
        {"metric": "single_copy_orthogroups", "value": len(files)},
        {"metric": "species", "value": len(species_order)},
        {"metric": "concatenated_alignment_aa", "value": len(next(iter(joined.values())))},
        {"metric": "species_names", "value": ",".join(species_order)},
        {"metric": "alignment_directory", "value": args.alignment_dir},
        {"metric": "output_fasta", "value": args.output_fasta},
    ]
    with args.summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["metric", "value"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
