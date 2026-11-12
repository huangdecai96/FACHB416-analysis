#!/usr/bin/env python3
"""Integrate miniprot and TBLASTN evidence for zero-member contracted families."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def interval_union(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = []
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def overlap_gene(
    chrom: str, start: int, end: int, genes: dict[str, list[dict[str, str]]]
) -> tuple[str, str, int, float, float]:
    hits = []
    locus_span = max(1, end - start + 1)
    for gene in genes.get(chrom, []):
        gs, ge = int(gene["start"]), int(gene["end"])
        overlap = max(0, min(end, ge) - max(start, gs) + 1)
        if overlap:
            gene_span = max(1, ge - gs + 1)
            hits.append(
                (
                    overlap,
                    overlap / locus_span,
                    overlap / gene_span,
                    gene["gene_id"],
                    gene["clean_protein_id"],
                )
            )
    if not hits:
        return "", "", 0, 0.0, 0.0
    overlap, locus_fraction, gene_fraction, gene_id, protein_id = max(hits)
    return gene_id, protein_id, overlap, locus_fraction, gene_fraction


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query-map", required=True)
    ap.add_argument("--family-table", required=True)
    ap.add_argument("--miniprot-paf", required=True)
    ap.add_argument("--tblastn", required=True)
    ap.add_argument("--gene-map", required=True)
    ap.add_argument("--orthogroups", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    qmeta = {r["query_id"]: r for r in read_tsv(Path(args.query_map))}
    families = {r["family_id"]: r for r in read_tsv(Path(args.family_table))}
    genes: dict[str, list[dict[str, str]]] = defaultdict(list)
    protein_to_og: dict[str, str] = {}
    for row in read_tsv(Path(args.gene_map)):
        genes[row["seqid"]].append(row)
    for row in read_tsv(Path(args.orthogroups)):
        for protein in [x.strip() for x in row.get("FCHB416", "").split(",") if x.strip()]:
            protein_to_og[protein] = row["Orthogroup"]

    miniprot_best: dict[str, dict[str, object]] = {}
    paf = Path(args.miniprot_paf)
    if paf.exists():
        with paf.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    continue
                qid, qlen, qstart, qend, strand, chrom, _, tstart, tend, matches, block, mapq = fields[:12]
                if qid not in qmeta:
                    continue
                tags = {x[:2]: x[5:] for x in fields[12:] if len(x) > 5 and x[2:5] in (":i:", ":f:", ":Z:")}
                score = float(tags.get("AS", matches))
                rec = {
                    "query_id": qid,
                    "chrom": chrom,
                    "start": int(tstart) + 1,
                    "end": int(tend),
                    "strand": strand,
                    "query_length": int(qlen),
                    "query_coverage": (int(qend) - int(qstart)) / int(qlen),
                    "identity": int(matches) / int(block) if int(block) else 0,
                    "mapq": int(mapq),
                    "score": score,
                    "frameshifts": int(tags.get("fs", 0)),
                    "premature_stops": int(tags.get("st", 0)),
                }
                if qid not in miniprot_best or score > float(miniprot_best[qid]["score"]):
                    miniprot_best[qid] = rec

    blast_by_query_locus: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    with Path(args.tblastn).open("r", encoding="utf-8", newline="") as handle:
        fields = ["query_id", "chrom", "pident", "length", "query_length", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        for row in csv.DictReader(handle, delimiter="\t", fieldnames=fields):
            if row["query_id"] in qmeta:
                strand = "+" if int(row["sstart"]) <= int(row["send"]) else "-"
                blast_by_query_locus[(row["query_id"], row["chrom"], strand)].append(row)
    blast_best: dict[str, dict[str, object]] = {}
    max_locus_gap = 50_000
    for (qid, chrom, strand), rows in blast_by_query_locus.items():
        ordered = sorted(rows, key=lambda r: min(int(r["sstart"]), int(r["send"])))
        clusters: list[list[dict[str, str]]] = []
        for row in ordered:
            start = min(int(row["sstart"]), int(row["send"]))
            if not clusters:
                clusters.append([row])
                continue
            prior_end = max(max(int(r["sstart"]), int(r["send"])) for r in clusters[-1])
            if start - prior_end <= max_locus_gap:
                clusters[-1].append(row)
            else:
                clusters.append([row])
        for locus_index, locus_rows in enumerate(clusters, 1):
            qlen = int(locus_rows[0]["query_length"])
            aligned = interval_union([(int(r["qstart"]), int(r["qend"])) for r in locus_rows])
            total_len = sum(int(r["length"]) for r in locus_rows)
            weighted_identity = sum(float(r["pident"]) * int(r["length"]) for r in locus_rows) / total_len
            rec = {
                "query_id": qid,
                "chrom": chrom,
                "start": min(min(int(r["sstart"]), int(r["send"])) for r in locus_rows),
                "end": max(max(int(r["sstart"]), int(r["send"])) for r in locus_rows),
                "strand": strand,
                "locus_index": locus_index,
                "query_coverage": min(1.0, aligned / qlen),
                "identity": weighted_identity / 100,
                "evalue": min(float(r["evalue"]) for r in locus_rows),
                "bitscore": sum(float(r["bitscore"]) for r in locus_rows),
            }
            rank = float(rec["query_coverage"]) * float(rec["identity"]) * float(rec["bitscore"])
            old = blast_best.get(qid)
            old_rank = (
                float(old["query_coverage"]) * float(old["identity"]) * float(old["bitscore"])
                if old
                else -1.0
            )
            if rank > old_rank:
                blast_best[qid] = rec

    query_rows = []
    family_evidence: dict[str, list[dict[str, object]]] = defaultdict(list)
    for qid, meta in qmeta.items():
        mp = miniprot_best.get(qid, {})
        bl = blast_best.get(qid, {})
        evidence = mp if float(mp.get("query_coverage", 0)) >= float(bl.get("query_coverage", 0)) else bl
        qcov = float(evidence.get("query_coverage", 0))
        ident = float(evidence.get("identity", 0))
        strong = qcov >= 0.70 and ident >= 0.25
        fragment = qcov >= 0.20 and ident >= 0.20
        intact_miniprot = (
            evidence is mp
            and bool(mp)
            and qcov >= 0.80
            and ident >= 0.30
            and int(mp.get("frameshifts", 0)) == 0
            and int(mp.get("premature_stops", 0)) == 0
        )
        chrom = str(evidence.get("chrom", ""))
        start, end = int(evidence.get("start", 0)), int(evidence.get("end", 0))
        gene_id, protein_id, overlap_bp, locus_fraction, gene_fraction = (
            overlap_gene(chrom, start, end, genes) if chrom else ("", "", 0, 0.0, 0.0)
        )
        meaningful_gene_overlap = overlap_bp >= 300 and locus_fraction >= 0.20 and gene_fraction >= 0.20
        if not meaningful_gene_overlap:
            gene_id, protein_id = "", ""
        row: dict[str, object] = {
            **meta,
            "best_method": "miniprot" if evidence is mp and mp else ("tblastn" if bl else "none"),
            "chrom": chrom,
            "start": start or "",
            "end": end or "",
            "query_coverage": qcov,
            "identity": ident,
            "strong_hit": "yes" if strong else "no",
            "fragment_hit": "yes" if fragment else "no",
            "intact_miniprot_hit": "yes" if intact_miniprot else "no",
            "mapq": evidence.get("mapq", ""),
            "alignment_score": evidence.get("score", evidence.get("bitscore", "")),
            "frameshifts": evidence.get("frameshifts", ""),
            "premature_stops": evidence.get("premature_stops", ""),
            "overlap_gene_id": gene_id,
            "overlap_clean_protein_id": protein_id,
            "overlap_gene_orthogroup": protein_to_og.get(protein_id, ""),
            "overlap_bp": overlap_bp,
            "locus_overlap_fraction": locus_fraction,
            "gene_overlap_fraction": gene_fraction,
        }
        query_rows.append(row)
        family_evidence[meta["family_id"]].append(row)

    family_rows = []
    for family, meta in families.items():
        evidence = family_evidence.get(family, [])
        strong = [r for r in evidence if r["strong_hit"] == "yes"]
        intact = [r for r in evidence if r["intact_miniprot_hit"] == "yes"]
        fragments = [r for r in evidence if r["fragment_hit"] == "yes"]
        annotated = [r for r in strong if r["overlap_gene_id"]]
        unannotated_intact = [r for r in intact if not r["overlap_gene_id"]]
        unannotated_other = [r for r in strong if not r["overlap_gene_id"] and r["intact_miniprot_hit"] != "yes"]
        species_support = len({str(r["source_species"]) for r in intact})
        best = max(evidence, key=lambda r: float(r["query_coverage"]) * float(r["identity"]), default={})
        if annotated:
            classification = "annotated_gene_or_orthogroup_assignment_rescue"
        elif unannotated_intact and (
            species_support >= 2
            or (float(best.get("query_coverage", 0)) >= 0.90 and float(best.get("identity", 0)) >= 0.40)
        ):
            classification = "unannotated_intact_locus_candidate"
        elif unannotated_intact or unannotated_other:
            classification = "unannotated_locus_low_support"
        elif fragments:
            classification = "fragment_or_pseudogene_candidate"
        else:
            classification = "no_detectable_genomic_homolog_candidate_absence"
        family_rows.append(
            {
                **meta,
                "rescue_class": classification,
                "queries_with_strong_hit": len(strong),
                "source_species_with_strong_hit": species_support,
                "queries_with_intact_miniprot_hit": len(intact),
                "strong_hits_overlapping_annotated_gene": len(annotated),
                "intact_hits_in_unannotated_region": len(unannotated_intact),
                "other_strong_hits_in_unannotated_region": len(unannotated_other),
                "queries_with_fragment_hit": len(fragments),
                "best_query_id": best.get("query_id", ""),
                "best_query_coverage": best.get("query_coverage", ""),
                "best_identity": best.get("identity", ""),
                "best_locus": f"{best.get('chrom', '')}:{best.get('start', '')}-{best.get('end', '')}" if best else "",
                "overlap_genes": ",".join(sorted({str(r["overlap_gene_id"]) for r in annotated if r["overlap_gene_id"]})),
                "overlap_gene_orthogroups": ",".join(sorted({str(r["overlap_gene_orthogroup"]) for r in annotated if r["overlap_gene_orthogroup"]})),
            }
        )

    for filename, rows in [("zero_member_query_evidence.tsv", query_rows), ("zero_member_family_rescue.tsv", family_rows)]:
        with (out / filename).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
