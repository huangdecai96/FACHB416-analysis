#!/usr/bin/env bash
set -euo pipefail

PROJ="/PROJECT_ROOT/Gene/paper2_finalref_rnaseq_20260713"
FINAL_FASTA="/PROJECT_ROOT/Gene/FCHB416_annotation/final_ref_qc_20260608/T_obliquus_FCHB416_chr17_Pt169551_Mt39134.fasta"
FINAL_GTF="/PROJECT_ROOT/Gene/FCHB416_annotation/braker3_out_FCHB416/braker.gtf"
SAMPLE_SHEET="$PROJ/metadata/paper2_finalref_sample_sheet_20260713.tsv"
THREADS="${THREADS:-48}"
SORT_THREADS="${SORT_THREADS:-8}"

mkdir -p "$PROJ"/{bam,counts,index,logs,metadata,qc,scripts,tmp}
export TMPDIR="$PROJ/tmp"

{
  date
  hostname
  df -h /PROJECT_ROOT
  echo "FINAL_FASTA=$FINAL_FASTA"
  echo "FINAL_GTF=$FINAL_GTF"
  echo "THREADS=$THREADS"
  echo "SORT_THREADS=$SORT_THREADS"
  command -v hisat2
  hisat2 --version | head -3
  command -v hisat2-build
  command -v samtools
  samtools --version | head -3
  command -v featureCounts
  featureCounts -v 2>&1 | head -3 || true
  command -v multiqc || true
} > "$PROJ/logs/paper2_finalref_recount_environment_20260713.log" 2>&1

test -s "$FINAL_FASTA"
test -s "$FINAL_GTF"
test -s "$SAMPLE_SHEET"

if [ ! -s "$FINAL_FASTA.fai" ]; then
  samtools faidx "$FINAL_FASTA"
fi
cp "$FINAL_FASTA.fai" "$PROJ/qc/final_reference.faidx.tsv"

if [ ! -s "$PROJ/index/FCHB416_finalref.1.ht2" ]; then
  hisat2-build -p "$THREADS" "$FINAL_FASTA" "$PROJ/index/FCHB416_finalref" \
    > "$PROJ/logs/hisat2_build_finalref.log" 2>&1
fi

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r sample condition concentration source read1 read2; do
  bam="$PROJ/bam/${sample}.finalref.sorted.bam"
  if [ -s "$bam" ] && [ -s "$bam.bai" ]; then
    echo "$(date) SKIP existing $sample"
  else
    echo "$(date) START $sample"
    test -s "$read1"
    test -s "$read2"
    hisat2 -p "$THREADS" --dta -x "$PROJ/index/FCHB416_finalref" -1 "$read1" -2 "$read2" \
      2> "$PROJ/logs/${sample}.hisat2.log" \
    | samtools sort -@ "$SORT_THREADS" -o "$bam" -
    samtools index -@ "$SORT_THREADS" "$bam"
    echo "$(date) DONE $sample"
  fi
  samtools flagstat -@ "$SORT_THREADS" "$bam" > "$PROJ/qc/${sample}.flagstat.txt"
  samtools idxstats "$bam" > "$PROJ/qc/${sample}.idxstats.tsv"
done | tee "$PROJ/logs/paper2_finalref_alignment_progress_20260713.log"

featureCounts -T "$THREADS" -s 0 -p --countReadPairs -t exon -g gene_id \
  -a "$FINAL_GTF" -o "$PROJ/counts/finalref_braker_strict_featureCounts.txt" "$PROJ"/bam/*.finalref.sorted.bam \
  > "$PROJ/logs/finalref_braker_strict_featureCounts.stdout.log" \
  2> "$PROJ/logs/finalref_braker_strict_featureCounts.stderr.log"

featureCounts -T "$THREADS" -s 0 -p --countReadPairs -M --fraction -t exon -g gene_id \
  -a "$FINAL_GTF" -o "$PROJ/counts/finalref_braker_multimap_fraction_featureCounts.txt" "$PROJ"/bam/*.finalref.sorted.bam \
  > "$PROJ/logs/finalref_braker_multimap_fraction_featureCounts.stdout.log" \
  2> "$PROJ/logs/finalref_braker_multimap_fraction_featureCounts.stderr.log"

multiqc "$PROJ/logs" "$PROJ/qc" "$PROJ/counts" -o "$PROJ/qc/multiqc_finalref" \
  > "$PROJ/logs/multiqc_finalref.stdout.log" 2> "$PROJ/logs/multiqc_finalref.stderr.log" || true

{
  date
  df -h /PROJECT_ROOT
  ls -lh "$PROJ/counts"
} > "$PROJ/logs/paper2_finalref_recount_finish_20260713.log" 2>&1

echo "paper2 finalref recount finished: $PROJ"
