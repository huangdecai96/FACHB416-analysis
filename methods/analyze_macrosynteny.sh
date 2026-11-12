#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-/PROJECT_ROOT/Gene/fachb416_integrated_synteny_allcomb_20260723}"
SRC="${SRC:-/PROJECT_ROOT/Gene/fachb416_integrated_synteny_20260723/prepared}"
THREADS="${THREADS:-32}"
MCS="${MCS:-/PROJECT_ROOT/miniconda3/envs/fchb416_v15_figures/bin/MCScanX}"

mkdir -p "$ROOT"/{logs,scripts,prepared,dbs,pairs,results,figures}
LOG="$ROOT/logs/run_integrated_multistrain_synteny_allcomb_20260723.log"
exec > >(tee -a "$LOG") 2>&1

echo "== start all-combination synteny =="
date -Is
echo "ROOT=$ROOT"
echo "SRC=$SRC"
echo "THREADS=$THREADS"
echo "MCS=$MCS"
echo "== user-space quota =="
quota -s || true
df -h /PROJECT_ROOT || true

SPECIES=(
  FCHB416
  Tobliquus_SAG276_10
  Tdistendus_CCAP276_37
  Araciborskii_CCAP276_35
  Ppectinatus_SAG2003
  Cmucronulatus_SAG2332
  Sobtusus_CCAP276_36
  Shindakii_SAG47_86
  Sproducto_capitatus_SAG21_81
)

ADJACENT_PAIRS=(
  FCHB416__Tobliquus_SAG276_10
  Tobliquus_SAG276_10__Tdistendus_CCAP276_37
  Tdistendus_CCAP276_37__Araciborskii_CCAP276_35
  Araciborskii_CCAP276_35__Ppectinatus_SAG2003
  Ppectinatus_SAG2003__Cmucronulatus_SAG2332
  Cmucronulatus_SAG2332__Sobtusus_CCAP276_36
  Sobtusus_CCAP276_36__Shindakii_SAG47_86
  Shindakii_SAG47_86__Sproducto_capitatus_SAG21_81
)

printf "%s\n" "${SPECIES[@]}" > "$ROOT/results/species_order_9way.txt"
printf "%s\n" "${ADJACENT_PAIRS[@]}" > "$ROOT/results/adjacent_pair_keys.txt"

echo "== copy prepared protein and gff files =="
for sp in "${SPECIES[@]}"; do
  test -s "$SRC/${sp}.pep.fa"
  test -s "$SRC/${sp}.gff"
  cp -f "$SRC/${sp}.pep.fa" "$ROOT/prepared/${sp}.pep.fa"
  cp -f "$SRC/${sp}.gff" "$ROOT/prepared/${sp}.gff"
  echo "$sp proteins=$(grep -c '^>' "$ROOT/prepared/${sp}.pep.fa") gff_lines=$(wc -l < "$ROOT/prepared/${sp}.gff")"
done

echo "== build DIAMOND databases once =="
for sp in "${SPECIES[@]}"; do
  if [ ! -s "$ROOT/dbs/${sp}.dmnd" ]; then
    diamond makedb --in "$ROOT/prepared/${sp}.pep.fa" -d "$ROOT/dbs/${sp}" --quiet
  fi
done

echo "== run all pairwise MCScanX combinations =="
for ((i=0; i<${#SPECIES[@]}; i++)); do
  for ((j=i+1; j<${#SPECIES[@]}; j++)); do
    A="${SPECIES[$i]}"
    B="${SPECIES[$j]}"
    PAIR="${A}_vs_${B}"
    PAIR_DIR="$ROOT/pairs/$PAIR"
    mkdir -p "$PAIR_DIR"
    echo "== pair $PAIR =="
    if [ -s "$PAIR_DIR/${PAIR}.collinearity" ]; then
      echo "skip existing complete collinearity: $PAIR_DIR/${PAIR}.collinearity"
      continue
    fi
    cat "$ROOT/prepared/${A}.gff" "$ROOT/prepared/${B}.gff" > "$PAIR_DIR/${PAIR}.gff"
    diamond blastp \
      --query "$ROOT/prepared/${A}.pep.fa" \
      --db "$ROOT/dbs/${B}.dmnd" \
      --out "$PAIR_DIR/${PAIR}.blast" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
      --evalue 1e-5 \
      --max-target-seqs 25 \
      --threads "$THREADS" \
      --quiet
    (
      cd "$PAIR_DIR"
      "$MCS" "$PAIR" > "MCScanX.stdout" 2> "MCScanX.stderr"
    )
  done
done

PY="$ROOT/scripts/parse_allcomb_mcscanx.py"
cat > "$PY" <<'PYCODE'
#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PREP = ROOT / "prepared"
PAIRS = ROOT / "pairs"


def read_list(path: Path) -> list[str]:
    return [x.strip() for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]


SPECIES = read_list(RESULTS / "species_order_9way.txt")
ADJ_KEYS = set(read_list(RESULTS / "adjacent_pair_keys.txt"))


def load_gff(species: str) -> dict[str, tuple[str, int, int]]:
    coords = {}
    with (PREP / f"{species}.gff").open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, gene, start, end = parts[:4]
            try:
                coords[gene] = (chrom, int(float(start)), int(float(end)))
            except ValueError:
                continue
    return coords


COORDS = {sp: load_gff(sp) for sp in SPECIES}


def parse_collinearity(path: Path, species_a: str, species_b: str) -> list[dict[str, object]]:
    rows = []
    block_id = ""
    score = ""
    evalue = ""
    npairs = ""
    header_re = re.compile(r"##\s+(Alignment\s+\d+):\s+score=([^e]+)e_value=([^\s]+)\s+N=([^\s]+)")
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("## Alignment"):
                m = header_re.search(line)
                if m:
                    block_id = m.group(1).strip()
                    score = m.group(2).strip()
                    evalue = m.group(3).strip()
                    npairs = m.group(4).strip()
                else:
                    block_id = line.lstrip("#").strip()
                continue
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            g1 = parts[1].strip()
            g2 = parts[2].strip()
            ca = COORDS[species_a].get(g1)
            cb = COORDS[species_b].get(g2)
            if ca is None or cb is None:
                ca2 = COORDS[species_a].get(g2)
                cb2 = COORDS[species_b].get(g1)
                if ca2 is not None and cb2 is not None:
                    g1, g2 = g2, g1
                    ca, cb = ca2, cb2
            if ca is None or cb is None:
                continue
            rows.append({
                "species_a": species_a,
                "species_b": species_b,
                "pair_key": f"{species_a}__{species_b}",
                "is_adjacent_chain_pair": f"{species_a}__{species_b}" in ADJ_KEYS,
                "block_id": block_id,
                "block_score": score,
                "block_evalue": evalue,
                "block_n_pairs": npairs,
                "gene_a": g1,
                "gene_b": g2,
                "chr_a": ca[0],
                "start_a": ca[1],
                "end_a": ca[2],
                "chr_b": cb[0],
                "start_b": cb[1],
                "end_b": cb[2],
            })
    return rows


all_rows: list[dict[str, object]] = []
summary_rows = []
chr_rows = []
for i, sp_a in enumerate(SPECIES):
    for sp_b in SPECIES[i+1:]:
        pair = f"{sp_a}_vs_{sp_b}"
        col = PAIRS / pair / f"{pair}.collinearity"
        status = "complete" if col.exists() and col.stat().st_size > 0 else "missing_collinearity"
        rows = parse_collinearity(col, sp_a, sp_b) if status == "complete" else []
        all_rows.extend(rows)
        blocks = {r["block_id"] for r in rows if r["block_id"]}
        summary_rows.append({
            "species_a": sp_a,
            "species_b": sp_b,
            "pair_key": f"{sp_a}__{sp_b}",
            "is_adjacent_chain_pair": f"{sp_a}__{sp_b}" in ADJ_KEYS,
            "mcscanx_status": status,
            "collinear_blocks": len(blocks),
            "collinear_gene_pairs": len(rows),
            "unique_genes_a": len({r["gene_a"] for r in rows}),
            "unique_genes_b": len({r["gene_b"] for r in rows}),
        })
        by_chr = defaultdict(int)
        for r in rows:
            by_chr[(r["chr_a"], r["chr_b"])] += 1
        for (ca, cb), count in sorted(by_chr.items()):
            chr_rows.append({
                "species_a": sp_a,
                "species_b": sp_b,
                "pair_key": f"{sp_a}__{sp_b}",
                "is_adjacent_chain_pair": f"{sp_a}__{sp_b}" in ADJ_KEYS,
                "chr_a": ca,
                "chr_b": cb,
                "collinear_gene_pairs": count,
            })


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


pair_fields = [
    "species_a", "species_b", "pair_key", "is_adjacent_chain_pair",
    "block_id", "block_score", "block_evalue", "block_n_pairs",
    "gene_a", "gene_b", "chr_a", "start_a", "end_a", "chr_b", "start_b", "end_b"
]
summary_fields = [
    "species_a", "species_b", "pair_key", "is_adjacent_chain_pair",
    "mcscanx_status", "collinear_blocks", "collinear_gene_pairs", "unique_genes_a", "unique_genes_b"
]
chr_fields = ["species_a", "species_b", "pair_key", "is_adjacent_chain_pair", "chr_a", "chr_b", "collinear_gene_pairs"]

write_tsv(RESULTS / "allcomb_collinear_gene_pairs_with_coords.tsv", all_rows, pair_fields)
write_tsv(RESULTS / "allcomb_mcscanx_summary.tsv", summary_rows, summary_fields)
write_tsv(RESULTS / "allcomb_chr_pair_matrix.tsv", chr_rows, chr_fields)

report = [
    "# 9-way close-relative all-combination MCScanX summary",
    "",
    f"- Species: {len(SPECIES)}",
    f"- Pair combinations: {len(summary_rows)}",
    f"- Adjacent chain pairs: {sum(1 for r in summary_rows if r['is_adjacent_chain_pair'])}",
    f"- Complete pairs: {sum(1 for r in summary_rows if r['mcscanx_status'] == 'complete')}",
    f"- Total collinear gene-pair rows: {len(all_rows)}",
    "",
    "This run extends the FACHB-centered synteny layer to all pairwise combinations among FACHB-416 and eight close relatives.",
]
(RESULTS / "allcomb_mcscanx_report_zh.md").write_text("\n".join(report) + "\n", encoding="utf-8")
PYCODE

python "$PY"

echo "== result summary =="
cat "$ROOT/results/allcomb_mcscanx_report_zh.md"
cat "$ROOT/results/allcomb_mcscanx_summary.tsv"
echo "== end all-combination synteny =="
date -Is
