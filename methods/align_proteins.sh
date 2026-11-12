#!/usr/bin/env bash
set -euo pipefail

# Code follows this file; ROOT continues to identify the external data project.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

USER_ROOT=/PROJECT_ROOT
ROOT=$USER_ROOT/Gene/paper1_round31_five_audits_20260719
SC=$ROOT/inputs/expanded20/orthofinder_expanded20_original/OrthoFinder/Results_Round31_original/Single_Copy_Orthologue_Sequences
OUT=$ROOT/results/expanded20_species_tree_ml_audit
ALIGN=$OUT/alignments
TRIM=$OUT/trimmed_alignments
PY=$USER_ROOT/miniconda3/envs/orthofinder_py310/bin/python
MAFFT=$USER_ROOT/miniconda3/envs/orthofinder_py310/bin/mafft
TRIMAL=$USER_ROOT/miniconda3/envs/orthofinder_py310/bin/trimal
IQTREE=$USER_ROOT/miniconda3/envs/orthofinder_py310/bin/iqtree3
RSCRIPT=$USER_ROOT/miniconda3/envs/R42/bin/Rscript
LOG=$ROOT/logs/22_expanded20_species_tree_ml_audit.log

mkdir -p "$ALIGN" "$TRIM"
exec > >(tee -a "$LOG") 2>&1
echo "run_start=$(date --iso-8601=seconds)"
echo "single_copy_fastas=$(find "$SC" -maxdepth 1 -type f -name '*.fa' | wc -l)"

export MAFFT TRIMAL ALIGN TRIM
find "$SC" -maxdepth 1 -type f -name '*.fa' -print0 | \
  xargs -0 -n 1 -P 24 bash -c '
    input="$1"
    stem=$(basename "$input" .fa)
    "$MAFFT" --auto --anysymbol --thread 1 --quiet "$input" > "$ALIGN/$stem.aln.fa"
    "$TRIMAL" -in "$ALIGN/$stem.aln.fa" -out "$TRIM/$stem.trim.fa" -automated1
  ' _

"$PY" "$SCRIPT_DIR/concatenate_single_copy_alignments.py" \
  --alignment-dir "$TRIM" \
  --output-fasta "$OUT/expanded20_490_singlecopy_concat.faa" \
  --partition-tsv "$OUT/expanded20_490_singlecopy_partitions.tsv" \
  --summary-tsv "$OUT/expanded20_490_singlecopy_alignment_summary.tsv"
cat "$OUT/expanded20_490_singlecopy_alignment_summary.tsv"

"$IQTREE" \
  -s "$OUT/expanded20_490_singlecopy_concat.faa" \
  -m MFP -mset LG,WAG,JTT,VT \
  -B 1000 --bnni -alrt 1000 \
  -T 24 --seed 416 \
  --prefix "$OUT/FCHB416_expanded20_490SC_ML"

if [[ ! -x "$RSCRIPT" ]]; then RSCRIPT=$(command -v Rscript); fi
"$RSCRIPT" "$SCRIPT_DIR/root_iqtree_species_tree.R" \
  "$OUT/FCHB416_expanded20_490SC_ML.treefile" \
  "$OUT/FCHB416_expanded20_490SC_ML.rooted.nwk"

touch "$OUT/expanded20_species_tree_ml_audit.complete"
echo "run_finish=$(date --iso-8601=seconds)"
quota -s || true
