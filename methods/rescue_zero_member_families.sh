#!/usr/bin/env bash
set -euo pipefail

# Code follows this file; ROOT continues to identify the external data project.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

USER_ROOT=/PROJECT_ROOT
OLD=$USER_ROOT/Gene/FCHB416_comparative_genomics_v4_20260611/comparative_genomics_nitrogen_v4
ROOT=$USER_ROOT/Gene/paper1_round31_five_audits_20260719
OUT=$ROOT/results/zero_member_rescue
LOG=$ROOT/logs/04_zero_member_rescue.log
PYTHON=$USER_ROOT/miniconda3/bin/python3
MINIPROT=$ROOT/envs/miniprot/bin/miniprot
BLAST_BIN=$USER_ROOT/miniconda3/envs/orthofinder_env/bin
GENOME=$OLD/01_inputs/FCHB416/FCHB416_nuclear17.unmasked.fa
OG=$OLD/01_inputs/orthofinder_10plus_species/OrthoFinder/Results_Jun12/Orthogroups/Orthogroups.tsv
CONTRACTED=$OLD/08_cafe5/cafe5_v12_timetree_verified_20260612/parsed/cafe_v12_FCHB41613_significant_contracted_families.tsv

mkdir -p "$OUT" "$OUT/blastdb"
exec > >(tee -a "$LOG") 2>&1
echo "run_start=$(date --iso-8601=seconds)"
echo "miniprot=$($MINIPROT --version 2>&1 || true)"
echo "blast=$($BLAST_BIN/tblastn -version 2>&1 | head -1)"
md5sum "$GENOME" "$OG" "$CONTRACTED"

"$PYTHON" "$SCRIPT_DIR/prepare_zero_member_queries.py" \
  --contracted "$CONTRACTED" \
  --orthogroups "$OG" \
  --proteome-dir "$OLD/01_inputs/orthofinder_10plus_species" \
  --out-dir "$OUT/queries"

"$BLAST_BIN/makeblastdb" -in "$GENOME" -dbtype nucl -parse_seqids -out "$OUT/blastdb/FCHB416_nuclear17" \
  > "$ROOT/logs/04_makeblastdb.log" 2>&1

"$MINIPROT" -I -t 32 "$GENOME" "$OUT/queries/zero_member_queries.faa" \
  > "$OUT/miniprot_zero_member.paf" 2> "$ROOT/logs/04_miniprot.log"

"$BLAST_BIN/tblastn" -query "$OUT/queries/zero_member_queries.faa" \
  -db "$OUT/blastdb/FCHB416_nuclear17" \
  -evalue 1e-5 -seg yes -soft_masking true -num_threads 32 -max_target_seqs 20 \
  -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' \
  -out "$OUT/tblastn_zero_member.tsv"

"$PYTHON" "$SCRIPT_DIR/classify_rescued_genes.py" \
  --query-map "$OUT/queries/zero_member_query_map.tsv" \
  --family-table "$OUT/queries/zero_member_families.tsv" \
  --miniprot-paf "$OUT/miniprot_zero_member.paf" \
  --tblastn "$OUT/tblastn_zero_member.tsv" \
  --gene-map "$OLD/01_inputs/FCHB416/FCHB416_gene_id_map.tsv" \
  --orthogroups "$OG" \
  --out-dir "$OUT/summary"

echo '[rescue_counts]'
cut -f6 "$OUT/summary/zero_member_family_rescue.tsv" | sort | uniq -c
echo "run_finish=$(date --iso-8601=seconds)"
