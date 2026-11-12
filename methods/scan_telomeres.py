#!/usr/bin/env python3
"""Scan scaffold ends for green-algal telomere repeats.

The script does not assume a fixed scaffold count.  It scans every FASTA
record and reports the first/last 20 kb telomere motif enrichment.
"""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path


MOTIFS = ("TTTAGGG", "CCCTAAA")


def read_fasta(path: Path):
    name = None
    chunks = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            yield name, "".join(chunks).upper()


def motif_stats(seq: str) -> tuple[int, int]:
    """Return total exact motif count and longest perfect tandem run."""
    total = sum(seq.count(motif) for motif in MOTIFS)
    longest = 0
    for motif in MOTIFS:
        pattern = re.compile(f"(?:{motif})+")
        for match in pattern.finditer(seq):
            longest = max(longest, len(match.group(0)) // len(motif))
    return total, longest


def classify_end(count: int, longest_run: int, min_count: int, min_run: int) -> bool:
    return count >= min_count or longest_run >= min_run


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--window", type=int, default=20_000)
    parser.add_argument("--major-cutoff", type=int, default=1_000_000)
    parser.add_argument("--min-count", type=int, default=5)
    parser.add_argument("--min-run", type=int, default=3)
    args = parser.parse_args()

    rows = []
    class_counts = Counter()
    major_rows = []

    for rank, (name, seq) in enumerate(read_fasta(args.fasta), start=1):
        length = len(seq)
        head = seq[: args.window]
        tail = seq[-args.window :] if length > args.window else seq
        head_count, head_run = motif_stats(head)
        tail_count, tail_run = motif_stats(tail)
        head_tel = classify_end(head_count, head_run, args.min_count, args.min_run)
        tail_tel = classify_end(tail_count, tail_run, args.min_count, args.min_run)

        if head_tel and tail_tel:
            tel_class = "double_telomere_complete_chromosome"
        elif head_tel or tail_tel:
            tel_class = "single_telomere_chromosome_arm"
        else:
            tel_class = "no_telomere"

        overlap_warning = "yes" if length < 2 * args.window else "no"
        row = {
            "rank": rank,
            "scaffold": name,
            "length_bp": length,
            "major_ge_1Mb": "yes" if length >= args.major_cutoff else "no",
            "head_telomere": "yes" if head_tel else "no",
            "tail_telomere": "yes" if tail_tel else "no",
            "head_motif_count": head_count,
            "head_longest_perfect_run": head_run,
            "tail_motif_count": tail_count,
            "tail_longest_perfect_run": tail_run,
            "end_windows_overlap": overlap_warning,
            "class": tel_class,
        }
        rows.append(row)
        class_counts[tel_class] += 1
        if length >= args.major_cutoff:
            major_rows.append(row)

    fields = [
        "rank",
        "scaffold",
        "length_bp",
        "major_ge_1Mb",
        "head_telomere",
        "tail_telomere",
        "head_motif_count",
        "head_longest_perfect_run",
        "tail_motif_count",
        "tail_longest_perfect_run",
        "end_windows_overlap",
        "class",
    ]
    with args.out_tsv.open("w", newline="") as out:
        out.write("\t".join(fields) + "\n")
        for row in rows:
            out.write("\t".join(str(row[field]) for field in fields) + "\n")

    with args.summary.open("w") as out:
        out.write("# Task 2 telomere motif scan summary\n")
        out.write(f"fasta\t{args.fasta}\n")
        out.write(f"scaffolds_total\t{len(rows)}\n")
        out.write(f"major_scaffolds_ge_1Mb\t{len(major_rows)}\n")
        out.write(f"window_bp\t{args.window}\n")
        out.write(f"motifs\t{','.join(MOTIFS)}\n")
        out.write(f"telomere_end_threshold\tcount>={args.min_count} OR longest_run>={args.min_run}\n\n")
        out.write("## Class counts, all scaffolds\n")
        for key in (
            "double_telomere_complete_chromosome",
            "single_telomere_chromosome_arm",
            "no_telomere",
        ):
            out.write(f"{key}\t{class_counts.get(key, 0)}\n")

        out.write("\n## Major scaffolds (>=1 Mb)\n")
        for row in major_rows:
            out.write(
                "{scaffold}\tlen={length_bp}\thead={head_telomere}"
                "(count={head_motif_count},run={head_longest_perfect_run})\t"
                "tail={tail_telomere}(count={tail_motif_count},run={tail_longest_perfect_run})\t"
                "{class}\n".format(**row)
            )

        out.write("\n## Telomere-positive scaffolds, all sizes\n")
        for row in rows:
            if row["head_telomere"] == "yes" or row["tail_telomere"] == "yes":
                out.write(
                    "{scaffold}\tlen={length_bp}\tmajor={major_ge_1Mb}\t"
                    "head={head_telomere}(count={head_motif_count},run={head_longest_perfect_run})\t"
                    "tail={tail_telomere}(count={tail_motif_count},run={tail_longest_perfect_run})\t"
                    "{class}\n".format(**row)
                )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
