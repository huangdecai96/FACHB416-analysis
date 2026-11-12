# File-interface adapter; statistical algorithm is unchanged.
aa <- grep('^--file=', commandArgs(FALSE), value=TRUE)
repo <- normalizePath(dirname(sub('^--file=', '', aa[1])), winslash='/')
args <- commandArgs(TRUE)
d <- if(length(args)>=1) args[1] else file.path(repo,'../FACHB416-source-data')
d <- normalizePath(d, winslash='/', mustWork=TRUE)
out <- if(length(args)>=2) args[2] else file.path(repo,'results/deseq2')
out <- normalizePath(out,winslash='/',mustWork=FALSE)
stopifnot(file.exists(file.path(d,'figure_data_map.tsv')))
stopifnot(out != d, !startsWith(out,paste0(d,'/')), !startsWith(d,paste0(out,'/')))
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop('Output folder is nonempty; use a new output folder')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(DESeq2))
input_dir<-file.path(d,'RNA');out_dir<-out
counts <- as.matrix(read.csv(file.path(input_dir, "paper2_finalref_braker_strict_raw_counts.csv"), row.names=1, check.names=FALSE))
storage.mode(counts) <- "numeric"
stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == round(counts)))
gene_info <- read.csv(file.path(input_dir, "paper2_finalref_gene_info.csv"), check.names=FALSE)
meta <- read.delim(file.path(input_dir, "sample_metadata.tsv"), check.names=FALSE)
rownames(meta) <- meta$sample
stopifnot(all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), , drop=FALSE]

meta$concentration <- factor(meta$concentration, levels = c("L", "H"))
meta$source <- factor(meta$source, levels = c("Nitrate", "Urea", "Ammonium"))
meta$condition <- factor(meta$condition, levels = c("LN", "LU", "LA", "HN", "HU", "HA"))

write.csv(counts, file.path(out_dir, "paper2_finalref_braker_strict_raw_counts.csv"))
write.csv(meta, file.path(out_dir, "paper2_finalref_metadata.csv"))
write.csv(gene_info, file.path(out_dir, "paper2_finalref_gene_info.csv"), row.names = FALSE)

keep <- rowSums(counts >= 10) >= 3
write.csv(
  data.frame(metric = c("genes_total", "genes_kept_count_ge10_in_ge3_samples"),
             value = c(nrow(counts), sum(keep))),
  file.path(out_dir, "paper2_finalref_filter_summary.csv"),
  row.names = FALSE
)

dds_factorial <- DESeqDataSetFromMatrix(
  countData = counts[keep, , drop = FALSE],
  colData = meta,
  design = ~ concentration + source + concentration:source
)
dds_factorial <- DESeq(dds_factorial)
saveRDS(dds_factorial, file.path(out_dir, "paper2_finalref_DESeq2_factorial_dds.rds"))
write.csv(data.frame(result_name = resultsNames(dds_factorial)),
          file.path(out_dir, "paper2_finalref_DESeq2_factorial_resultsNames.csv"),
          row.names = FALSE)

for (nm in resultsNames(dds_factorial)[-1]) {
  res <- results(dds_factorial, name = nm, alpha = 0.05)
  out <- as.data.frame(res)
  out$gene_id <- rownames(out)
  out <- out[, c("gene_id", setdiff(colnames(out), "gene_id"))]
  write.csv(out, file.path(out_dir, paste0("paper2_finalref_factorial_", make.names(nm), ".csv")),
            row.names = FALSE)
}

dds_lrt <- DESeq(
  dds_factorial,
  test = "LRT",
  reduced = ~ concentration + source
)
res_lrt <- as.data.frame(results(dds_lrt, alpha = 0.05))
res_lrt$gene_id <- rownames(res_lrt)
res_lrt <- res_lrt[, c("gene_id", setdiff(colnames(res_lrt), "gene_id"))]
write.csv(res_lrt, file.path(out_dir, "paper2_finalref_interaction_LRT_full_vs_no_interaction.csv"),
          row.names = FALSE)

dds_condition <- DESeqDataSetFromMatrix(
  countData = counts[keep, , drop = FALSE],
  colData = meta,
  design = ~ condition
)
dds_condition <- DESeq(dds_condition)
saveRDS(dds_condition, file.path(out_dir, "paper2_finalref_DESeq2_condition_dds.rds"))
write.csv(data.frame(result_name = resultsNames(dds_condition)),
          file.path(out_dir, "paper2_finalref_DESeq2_condition_resultsNames.csv"),
          row.names = FALSE)

contrast_list <- list(
  HN_vs_LN = c("condition", "HN", "LN"),
  HU_vs_LU = c("condition", "HU", "LU"),
  HA_vs_LA = c("condition", "HA", "LA"),
  HU_vs_HN = c("condition", "HU", "HN"),
  HA_vs_HN = c("condition", "HA", "HN"),
  HA_vs_HU = c("condition", "HA", "HU")
)

summary_rows <- list()
for (nm in names(contrast_list)) {
  res <- results(dds_condition, contrast = contrast_list[[nm]], alpha = 0.05)
  out <- as.data.frame(res)
  out$gene_id <- rownames(out)
  out <- out[, c("gene_id", setdiff(colnames(out), "gene_id"))]
  out <- out[order(out$padj, na.last = TRUE), ]
  write.csv(out, file.path(out_dir, paste0("paper2_finalref_condition_", nm, ".csv")),
            row.names = FALSE)
  sig <- subset(out, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
  write.csv(sig, file.path(out_dir, paste0("paper2_finalref_condition_", nm, "_sig_padj0.05_lfc1.csv")),
            row.names = FALSE)
  summary_rows[[nm]] <- data.frame(
    contrast = nm,
    tested_genes = nrow(out),
    sig_genes = nrow(sig),
    up_genes = sum(sig$log2FoldChange > 0),
    down_genes = sum(sig$log2FoldChange < 0)
  )
}
write.csv(do.call(rbind, summary_rows),
          file.path(out_dir, "paper2_finalref_condition_DEG_summary_padj0.05_lfc1.csv"),
          row.names = FALSE)

vsd <- vst(dds_condition, blind = FALSE)
vmat <- assay(vsd)
write.csv(vmat, file.path(out_dir, "paper2_finalref_vst_matrix.csv"))
write.csv(cor(vmat), file.path(out_dir, "paper2_finalref_vst_sample_correlation.csv"))

pca <- prcomp(t(vmat), scale. = FALSE)
pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 3)
pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC1_percent = pct[1],
  PC2_percent = pct[2],
  meta[rownames(pca$x), c("condition", "concentration", "source"), drop = FALSE]
)
write.csv(pca_df, file.path(out_dir, "paper2_finalref_vst_PCA_PC1_PC2.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "paper2_finalref_DESeq2_sessionInfo.txt"))
