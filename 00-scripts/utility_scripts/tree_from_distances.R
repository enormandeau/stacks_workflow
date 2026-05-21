#!/usr/bin/env Rscript

# Build phylogenetic tree from wide distances format

# Load libraries
library(ape)
library(phangorn)

# Parse user input
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = paste0(input_file, "_tree.png")

# Load data
distances = read.table(input_file, header=T, row.names=1, stringsAsFactors=F)
distances = as.matrix(distances)
distances2 = as.dist(distances)

# Build tree
tree = ape::nj(distances)
#tree = phangorn::upgma(distances)

png(output_file, width=800, height=20 * nrow(distances))
    plot(tree, main=input_file)
dev.off()
