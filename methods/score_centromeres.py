#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
from collections import Counter, defaultdict


def read_fai(path):
    lengths = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def write_tsv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_gff_attrs(attr):
    out = {}
    for item in attr.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
        elif " " in item:
            k, v = item.split(" ", 1)
            v = v.strip('"')
        else:
            continue
        out[k] = v
    return out


def parse_barrnap(path):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_gff_attrs(parts[8])
            product = attrs.get("product", attrs.get("Name", "rRNA"))
            rows.append({
                "source": "barrnap",
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "type": "rRNA",
                "subtype": product,
                "score": parts[5],
                "evalue": "",
                "note": "barrnap eukaryotic rRNA prediction",
            })
    return rows


def parse_trnascan(path):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("Sequence") or line.startswith("Name") or line.startswith("-"):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue
            strand = "+"
            if start > end:
                start, end = end, start
                strand = "-"
            aa_type = parts[4]
            codon = parts[5]
            score = parts[8]
            note = ""
            if len(parts) > 9:
                note = " ".join(parts[9:])
            rows.append({
                "source": "tRNAscan-SE",
                "seqid": parts[0],
                "start": start,
                "end": end,
                "strand": strand,
                "type": "tRNA",
                "subtype": f"tRNA-{aa_type}({codon})",
                "score": score,
                "evalue": "",
                "note": note,
            })
    return rows


def infer_rfam_class(model_name, description):
    text = f"{model_name} {description}".lower()
    if "rrna" in text or "ribosomal rna" in text or "5s" in text or "5.8s" in text or "18s" in text or "28s" in text:
        return "rRNA"
    if "trna" in text or "transfer rna" in text:
        return "tRNA"
    if "snorna" in text or "snoRNA".lower() in text or "small nucleolar" in text:
        return "snoRNA"
    if "snrna" in text or "spliceosomal" in text or "u1" in text or "u2" in text or "u4" in text or "u5" in text or "u6" in text:
        return "snRNA"
    if "mir" in text or "microrna" in text or "micro-rna" in text:
        return "miRNA-like"
    if "rnase p" in text:
        return "RNaseP_RNA"
    if "srp" in text or "signal recognition particle" in text:
        return "SRP_RNA"
    if "lncrna" in text or "long non-coding" in text:
        return "lncRNA-like"
    return "other_Rfam_ncRNA"


def parse_rfam_tblout(path, seq_lengths):
    if not os.path.exists(path):
        return []
    rows = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split(maxsplit=27)
            if len(parts) < 18:
                continue
            # Infernal cmscan --fmt 2 tblout columns:
            # idx, target name, accession, query name, accession, clan name,
            # mdl, mdl from, mdl to, seq from, seq to, strand, ...
            # Older/other formats may omit idx; support both layouts.
            has_idx = parts[0].isdigit()
            if has_idx:
                if len(parts) < 19:
                    continue
                model = parts[1]
                accession = parts[2]
                seqid = parts[3]
                seq_from_i, seq_to_i = 9, 10
                strand_i, score_i, evalue_i, inc_i = 11, 16, 17, 18
                desc_i = 27
            else:
                model = parts[0]
                accession = parts[1]
                seqid = parts[2]
                seq_from_i, seq_to_i = 7, 8
                strand_i, score_i, evalue_i, inc_i = 9, 14, 15, 16
                desc_i = 17
            if seqid not in seq_lengths:
                continue
            try:
                start = int(parts[seq_from_i])
                end = int(parts[seq_to_i])
            except (ValueError, IndexError):
                continue
            strand = parts[strand_i] if len(parts) > strand_i else ""
            score = parts[score_i] if len(parts) > score_i else ""
            evalue = parts[evalue_i] if len(parts) > evalue_i else ""
            inc = parts[inc_i] if len(parts) > inc_i else ""
            description = parts[desc_i] if len(parts) > desc_i else ""
            if start > end:
                start, end = end, start
            rows.append({
                "source": "Rfam/cmscan",
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand,
                "type": infer_rfam_class(model, description),
                "subtype": model,
                "score": score,
                "evalue": evalue,
                "note": f"{accession}; inc={inc}; {description}",
            })
    return rows


def parse_gene_gff3(path):
    genes = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue
            attrs = parse_gff_attrs(parts[8])
            genes.append({
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "id": attrs.get("ID", attrs.get("gene_id", "")),
            })
    return genes


def parse_edta_gff3(path):
    intervals = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_gff_attrs(parts[8])
            cls = attrs.get("Classification", attrs.get("Name", attrs.get("ID", parts[2])))
            intervals.append({
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "class": cls,
            })
    return intervals


def parse_trf_dat(path):
    intervals = []
    current = None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("Sequence:"):
                current = line.split(":", 1)[1].strip().split()[0]
                continue
            if current is None or not line[0].isdigit():
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            try:
                start = int(parts[0])
                end = int(parts[1])
                period = int(float(parts[2]))
                copies = float(parts[3])
                consensus = parts[13]
            except ValueError:
                continue
            intervals.append({
                "seqid": current,
                "start": start,
                "end": end,
                "period": period,
                "copies": copies,
                "consensus": consensus,
            })
    return intervals


def merge_intervals(intervals):
    by_seq = defaultdict(list)
    for rec in intervals:
        if rec["end"] < rec["start"]:
            rec["start"], rec["end"] = rec["end"], rec["start"]
        by_seq[rec["seqid"]].append((rec["start"], rec["end"]))
    merged = {}
    for seqid, vals in by_seq.items():
        vals.sort()
        out = []
        for s, e in vals:
            if not out or s > out[-1][1] + 1:
                out.append([s, e])
            else:
                out[-1][1] = max(out[-1][1], e)
        merged[seqid] = [(s, e) for s, e in out]
    return merged


def overlap_bp(intervals, start, end):
    total = 0
    for s, e in intervals:
        if e < start:
            continue
        if s > end:
            break
        total += max(0, min(e, end) - max(s, start) + 1)
    return total


def quantile(vals, q):
    vals = sorted(v for v in vals if not math.isnan(v))
    if not vals:
        return float("nan")
    pos = (len(vals) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[lo]
    return vals[lo] * (hi - pos) + vals[hi] * (pos - lo)


def minmax_scale(value, vals, invert=False):
    vals = [v for v in vals if not math.isnan(v)]
    if not vals:
        return 0.0
    mn = min(vals)
    mx = max(vals)
    if mx == mn:
        scaled = 0.0
    else:
        scaled = (value - mn) / (mx - mn)
    return 1.0 - scaled if invert else scaled


def build_windows(lengths, genes, edta, trf, softmask_seqs, window_size=250000):
    merged_te = merge_intervals([{"seqid": r["seqid"], "start": r["start"], "end": r["end"]} for r in edta])
    merged_trf = merge_intervals([{"seqid": r["seqid"], "start": r["start"], "end": r["end"]} for r in trf])
    genes_by_seq = defaultdict(list)
    for g in genes:
        genes_by_seq[g["seqid"]].append((g["start"], g["end"]))
    rows = []
    for seqid, length in lengths.items():
        seq = softmask_seqs.get(seqid, "")
        for start in range(1, length + 1, window_size):
            end = min(length, start + window_size - 1)
            span = end - start + 1
            gene_count = sum(1 for s, e in genes_by_seq.get(seqid, []) if not (e < start or s > end))
            te_bp = overlap_bp(merged_te.get(seqid, []), start, end)
            tandem_bp = overlap_bp(merged_trf.get(seqid, []), start, end)
            if seq:
                frag = seq[start - 1:end]
                soft_bp = sum(1 for ch in frag if ch.islower())
            else:
                soft_bp = 0
            internal = start > 500000 and end < (length - 500000)
            rows.append({
                "seqid": seqid,
                "start": start,
                "end": end,
                "midpoint_mb": (start + end) / 2 / 1_000_000,
                "window_bp": span,
                "gene_count": gene_count,
                "gene_density_per_mb": gene_count / span * 1_000_000,
                "te_bp": te_bp,
                "te_percent": te_bp / span * 100,
                "tandem_bp": tandem_bp,
                "tandem_percent": tandem_bp / span * 100,
                "softmask_bp": soft_bp,
                "softmask_percent": soft_bp / span * 100,
                "internal_window": "yes" if internal else "no",
            })
    return rows


def call_candidate_regions(windows):
    internal = [w for w in windows if w["internal_window"] == "yes"]
    all_gene = [w["gene_density_per_mb"] for w in internal]
    all_te = [w["te_percent"] for w in internal]
    all_tandem = [w["tandem_percent"] for w in internal]
    all_soft = [w["softmask_percent"] for w in internal]

    for w in windows:
        if w["internal_window"] != "yes":
            w["centromere_like_score"] = float("nan")
            continue
        score = (
            0.35 * minmax_scale(w["tandem_percent"], all_tandem) +
            0.25 * minmax_scale(w["te_percent"], all_te) +
            0.25 * minmax_scale(w["softmask_percent"], all_soft) +
            0.15 * minmax_scale(w["gene_density_per_mb"], all_gene, invert=True)
        )
        w["centromere_like_score"] = score

    by_chr = defaultdict(list)
    for w in windows:
        if w["internal_window"] == "yes":
            by_chr[w["seqid"]].append(w)

    candidates = []
    for seqid, vals in by_chr.items():
        vals.sort(key=lambda x: x["centromere_like_score"], reverse=True)
        if not vals:
            continue
        top = vals[0]
        support = []
        chr_tandem_q90 = quantile([v["tandem_percent"] for v in vals], 0.90)
        chr_te_q90 = quantile([v["te_percent"] for v in vals], 0.90)
        chr_soft_q90 = quantile([v["softmask_percent"] for v in vals], 0.90)
        chr_gene_q25 = quantile([v["gene_density_per_mb"] for v in vals], 0.25)
        if top["tandem_percent"] >= chr_tandem_q90:
            support.append("chr_top10_tandem")
        if top["te_percent"] >= chr_te_q90:
            support.append("chr_top10_TE")
        if top["softmask_percent"] >= chr_soft_q90:
            support.append("chr_top10_softmask")
        if top["gene_density_per_mb"] <= chr_gene_q25:
            support.append("chr_low25_gene_density")
        if top["tandem_percent"] >= 5:
            support.append("tandem_ge_5pct")
        confidence = "moderate"
        if len(support) >= 4:
            confidence = "higher"
        elif len(support) <= 1:
            confidence = "low"
        candidates.append({
            "seqid": seqid,
            "start": top["start"],
            "end": top["end"],
            "midpoint_mb": round(top["midpoint_mb"], 4),
            "centromere_like_score": round(top["centromere_like_score"], 4),
            "gene_density_per_mb": round(top["gene_density_per_mb"], 4),
            "te_percent": round(top["te_percent"], 4),
            "tandem_percent": round(top["tandem_percent"], 4),
            "softmask_percent": round(top["softmask_percent"], 4),
            "confidence": confidence,
            "support_flags": ",".join(support),
            "interpretation": "Exploratory centromere-like candidate; requires centromere-specific validation",
        })
    candidates.sort(key=lambda x: x["seqid"])
    return candidates


def trf_consensus_summary(trf_rows, candidates):
    out = []
    for cand in candidates:
        hits = [
            r for r in trf_rows
            if r["seqid"] == cand["seqid"] and not (r["end"] < cand["start"] or r["start"] > cand["end"])
        ]
        by_consensus = Counter()
        by_period = Counter()
        for h in hits:
            bp = max(0, min(h["end"], cand["end"]) - max(h["start"], cand["start"]) + 1)
            if bp <= 0:
                continue
            by_consensus[h["consensus"]] += bp
            by_period[h["period"]] += bp
        for consensus, bp in by_consensus.most_common(5):
            out.append({
                "seqid": cand["seqid"],
                "candidate_start": cand["start"],
                "candidate_end": cand["end"],
                "repeat_consensus": consensus,
                "repeat_bp_in_candidate": bp,
                "repeat_period_bp": "",
            })
        for period, bp in by_period.most_common(3):
            out.append({
                "seqid": cand["seqid"],
                "candidate_start": cand["start"],
                "candidate_end": cand["end"],
                "repeat_consensus": f"period_{period}",
                "repeat_bp_in_candidate": bp,
                "repeat_period_bp": period,
            })
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--softmask-fasta", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--trnascan", required=True)
    parser.add_argument("--barrnap", required=True)
    parser.add_argument("--rfam-tblout", required=True)
    parser.add_argument("--gene-gff3", required=True)
    parser.add_argument("--edta-gff3", required=True)
    parser.add_argument("--trf-dat", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    lengths = read_fai(args.fai)
    softmask_seqs = read_fasta(args.softmask_fasta)
    nc_rows = []
    nc_rows.extend(parse_trnascan(args.trnascan))
    nc_rows.extend(parse_barrnap(args.barrnap))
    nc_rows.extend(parse_rfam_tblout(args.rfam_tblout, lengths))

    nc_fields = ["source", "seqid", "start", "end", "strand", "type", "subtype", "score", "evalue", "note"]
    write_tsv(os.path.join(args.outdir, "FCHB416_nuclear_ncRNA_annotations.tsv"), nc_rows, nc_fields)

    summary = []
    for key, count in sorted(Counter((r["source"], r["type"]) for r in nc_rows).items()):
        summary.append({"source": key[0], "type": key[1], "count": count})
    write_tsv(os.path.join(args.outdir, "FCHB416_nuclear_ncRNA_summary_by_source_type.tsv"), summary, ["source", "type", "count"])

    chr_summary = []
    by_chr_type = Counter((r["seqid"], r["type"]) for r in nc_rows)
    for seqid in lengths:
        for nc_type in sorted(set(r["type"] for r in nc_rows)):
            chr_summary.append({"seqid": seqid, "type": nc_type, "count": by_chr_type[(seqid, nc_type)]})
    write_tsv(os.path.join(args.outdir, "FCHB416_nuclear_ncRNA_counts_by_chromosome.tsv"), chr_summary, ["seqid", "type", "count"])

    genes = parse_gene_gff3(args.gene_gff3)
    edta = parse_edta_gff3(args.edta_gff3)
    trf = parse_trf_dat(args.trf_dat)
    windows = build_windows(lengths, genes, edta, trf, softmask_seqs)
    candidates = call_candidate_regions(windows)
    consensus = trf_consensus_summary(trf, candidates)

    window_fields = [
        "seqid", "start", "end", "midpoint_mb", "window_bp", "gene_count", "gene_density_per_mb",
        "te_bp", "te_percent", "tandem_bp", "tandem_percent", "softmask_bp", "softmask_percent",
        "internal_window", "centromere_like_score"
    ]
    write_tsv(os.path.join(args.outdir, "FCHB416_centromere_like_window_evidence_250kb.tsv"), windows, window_fields)
    candidate_fields = [
        "seqid", "start", "end", "midpoint_mb", "centromere_like_score", "gene_density_per_mb",
        "te_percent", "tandem_percent", "softmask_percent", "confidence", "support_flags", "interpretation"
    ]
    write_tsv(os.path.join(args.outdir, "FCHB416_candidate_centromere_like_regions.tsv"), candidates, candidate_fields)
    write_tsv(
        os.path.join(args.outdir, "FCHB416_candidate_centromere_TRF_consensus_summary.tsv"),
        consensus,
        ["seqid", "candidate_start", "candidate_end", "repeat_consensus", "repeat_bp_in_candidate", "repeat_period_bp"],
    )

    te_class_counts = Counter()
    te_class_bp = Counter()
    for r in edta:
        cls = r["class"] or "unknown"
        major = cls.split("/")[0].split(":")[0]
        te_class_counts[major] += 1
        te_class_bp[major] += r["end"] - r["start"] + 1
    te_rows = [
        {"te_major_class": k, "feature_count": te_class_counts[k], "bp_sum_unmerged": te_class_bp[k]}
        for k in sorted(te_class_counts)
    ]
    write_tsv(os.path.join(args.outdir, "FCHB416_EDTA_TE_major_class_summary.tsv"), te_rows, ["te_major_class", "feature_count", "bp_sum_unmerged"])

    overall = [
        {"metric": "nuclear_chromosomes", "value": len(lengths)},
        {"metric": "nuclear_length_bp", "value": sum(lengths.values())},
        {"metric": "tRNAscan_tRNA_count", "value": sum(1 for r in nc_rows if r["source"] == "tRNAscan-SE")},
        {"metric": "barrnap_rRNA_count", "value": sum(1 for r in nc_rows if r["source"] == "barrnap")},
        {"metric": "Rfam_cmscan_count", "value": sum(1 for r in nc_rows if r["source"] == "Rfam/cmscan")},
        {"metric": "TRF_tandem_repeat_features", "value": len(trf)},
        {"metric": "EDTA_TE_features", "value": len(edta)},
        {"metric": "candidate_centromere_like_regions", "value": len(candidates)},
    ]
    write_tsv(os.path.join(args.outdir, "FCHB416_ncRNA_centromere_overall_summary.tsv"), overall, ["metric", "value"])


if __name__ == "__main__":
    main()
