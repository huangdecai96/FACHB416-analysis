#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: build_expanded20_time_trees.R <species_tree> <calibrations.tsv> <out_dir>", call. = FALSE)
}

tree_file <- args[[1]]
calibration_file <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

alias_map <- c(
  "Chlamydomonas reinhardtii" = "Creinhardtii",
  "Volvox carteri" = "Vcarteri",
  "Chlorella variabilis" = "Cvariabilis",
  "Ostreococcus tauri" = "Otauri",
  "Micromonas" = "Micromonas_RCC299",
  "Micromonas RCC299" = "Micromonas_RCC299"
)

resolve_tip <- function(x, tips) {
  x <- trimws(x)
  candidates <- unique(c(x, unname(alias_map[x])))
  candidates <- candidates[!is.na(candidates)]
  hit <- candidates[candidates %in% tips]
  if (length(hit) > 0) return(hit[[1]])
  normalized <- gsub("[^A-Za-z0-9]", "", tolower(x))
  normalized_tips <- gsub("[^A-Za-z0-9]", "", tolower(tips))
  idx <- match(normalized, normalized_tips)
  if (is.na(idx)) return(NA_character_)
  tips[[idx]]
}

tree <- read.tree(tree_file)
if (is.null(tree$edge.length)) stop("Input tree has no branch lengths", call. = FALSE)
outgroups <- c("Micromonas_RCC299", "Otauri")
missing_outgroups <- setdiff(outgroups, tree$tip.label)
if (length(missing_outgroups) > 0) {
  stop(paste("Missing outgroups:", paste(missing_outgroups, collapse = ",")), call. = FALSE)
}
tree <- root(tree, outgroup = outgroups, resolve.root = TRUE)
tree <- multi2di(tree, random = FALSE)
tree$node.label <- NULL
bad_edges <- which(is.na(tree$edge.length) | tree$edge.length <= 0)
if (length(bad_edges) > 0) tree$edge.length[bad_edges] <- 1e-8
write.tree(tree, file.path(out_dir, "expanded20_rooted_binary_input.nwk"))

cal <- read.delim(calibration_file, stringsAsFactors = FALSE, check.names = FALSE)
cal <- cal[cal$use_for_initial_time_tree %in% c("yes", "YES", "Yes", TRUE), , drop = FALSE]
if (nrow(cal) == 0) stop("No enabled calibrations", call. = FALSE)

cal_rows <- vector("list", nrow(cal))
for (i in seq_len(nrow(cal))) {
  tip_a <- resolve_tip(cal$taxon_a[[i]], tree$tip.label)
  tip_b <- resolve_tip(cal$taxon_b[[i]], tree$tip.label)
  if (is.na(tip_a) || is.na(tip_b)) {
    stop(paste("Unresolved calibration:", cal$calibration_id[[i]]), call. = FALSE)
  }
  node <- getMRCA(tree, c(tip_a, tip_b))
  cal_rows[[i]] <- data.frame(
    calibration_id = cal$calibration_id[[i]],
    taxon_a_tip = tip_a,
    taxon_b_tip = tip_b,
    node = node,
    age.min = as.numeric(cal$min_mya[[i]]),
    age.max = as.numeric(cal$max_mya[[i]]),
    source_status = cal$source_status[[i]],
    stringsAsFactors = FALSE
  )
}
cal_applied <- do.call(rbind, cal_rows)
write_tsv(cal_applied, file.path(out_dir, "expanded20_calibrations_applied.tsv"))
chronos_cal <- cal_applied[, c("node", "age.min", "age.max")]

variants <- data.frame(
  variant = c("correlated_l1", "correlated_l01", "relaxed_l1"),
  model = c("correlated", "correlated", "relaxed"),
  lambda = c(1, 0.1, 1),
  stringsAsFactors = FALSE
)

summary_rows <- list()
check_rows <- list()
for (i in seq_len(nrow(variants))) {
  v <- variants[i, ]
  fit <- chronos(
    tree,
    lambda = v$lambda,
    model = v$model,
    calibration = chronos_cal,
    control = chronos.control(iter.max = 2000, eval.max = 2000)
  )
  if (!is.ultrametric(fit, tol = 1e-6)) {
    stop(paste("Non-ultrametric output for", v$variant), call. = FALSE)
  }
  if (any(!is.finite(fit$edge.length)) || any(fit$edge.length < 0)) {
    stop(paste("Invalid branches for", v$variant), call. = FALSE)
  }
  out_tree <- file.path(out_dir, paste0("expanded20_", v$variant, ".nwk"))
  write.tree(fit, out_tree)
  depth <- node.depth.edgelength(fit)
  root_age <- max(depth[seq_along(fit$tip.label)])
  zero_edges <- sum(fit$edge.length == 0)
  summary_rows[[i]] <- data.frame(
    variant = v$variant,
    model = v$model,
    lambda = v$lambda,
    tips = length(fit$tip.label),
    rooted = is.rooted(fit),
    binary = is.binary(fit),
    ultrametric = is.ultrametric(fit, tol = 1e-6),
    root_age_mya = root_age,
    min_branch = min(fit$edge.length),
    zero_edges = zero_edges,
    tree_file = out_tree,
    stringsAsFactors = FALSE
  )
  for (j in seq_len(nrow(cal_applied))) {
    inferred <- root_age - depth[cal_applied$node[[j]]]
    check_rows[[length(check_rows) + 1]] <- data.frame(
      variant = v$variant,
      calibration_id = cal_applied$calibration_id[[j]],
      inferred_age_mya = inferred,
      min_mya = cal_applied$age.min[[j]],
      max_mya = cal_applied$age.max[[j]],
      within_range = inferred >= cal_applied$age.min[[j]] && inferred <= cal_applied$age.max[[j]],
      stringsAsFactors = FALSE
    )
  }
}

write_tsv(do.call(rbind, summary_rows), file.path(out_dir, "expanded20_time_tree_summary.tsv"))
write_tsv(do.call(rbind, check_rows), file.path(out_dir, "expanded20_calibration_checks.tsv"))
sink(file.path(out_dir, "sessionInfo.txt"))
sessionInfo()
sink()

