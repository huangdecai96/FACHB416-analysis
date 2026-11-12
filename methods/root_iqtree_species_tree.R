suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: root_iqtree_species_tree.R input.treefile output.rooted.nwk")
}

tree <- read.tree(args[[1]])
outgroup <- c("Otauri", "Micromonas_RCC299")
if (!all(outgroup %in% tree$tip.label)) {
  stop("Required outgroup taxa are absent")
}
rooted <- root(tree, outgroup = outgroup, resolve.root = TRUE)
write.tree(rooted, file = args[[2]])
