# Upstream scientific methods

These are the selected scientific methods behind the final data, retained for provenance and re-use with the required external inputs and tools. They are not standalone replacements for the portable replay scripts in the repository root.

The alignment, gene prediction evidence, centromere/telomere, phylogenetic calibration, synteny, localization and annotation-rescue methods are retained. The CAFE subdirectory contains model runs, result parsing and the final consensus summarizer. Original upstream project-relative input paths are retained where needed, with local machine roots replaced by /PROJECT_ROOT/.

The old server packaging driver, duplicate CAFE summarizer, duplicate DESeq2 implementation, historical figure wrappers and verification-result copies are excluded from this upload version. Run reproduce_deseq2.R for the verified final-count differential-expression workflow.

The RNA alignment wrapper is supplied once as align_rna.sh: its reference includes the nuclear chromosomes and both organelles. The protein-alignment and rescue wrappers include their actual query preparation, concatenation and rooting helpers. Their shell calls now locate bundled helpers beside the entry script; the scientific arguments and thresholds are unchanged. The MCScanX parser is already embedded in analyze_macrosynteny.sh and is not duplicated.

CAFE entry scripts take the external data-project directory as the first argument. That project must provide inputs/cafe_counts/{original,augmented}_counts.tsv and inputs/trees/{correlated_l1,correlated_l01,relaxed_l1}.nwk. Matching tables and trees are in the separate source-data archive under Comparative_genomics/. External binaries, upstream proteomes, raw FASTQ files and preprocessed alignment inputs must be supplied for full upstream reruns. These upstream workflows were source-checked, not recomputed during this packaging task.
