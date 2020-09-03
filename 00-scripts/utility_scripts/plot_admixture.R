#!/usr/bin/env Rscript

# Parse user input
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = paste0(input_file, ".png")

# Load libraries
library(adegenet)

# Load data
data <- read.table(input_file)
row.names(data) <- 1:nrow(data)
Qmatrix <- as.matrix(data)

# Plot figure
png(output_file, width=1800, height=300)
    compoplot(Qmatrix, legend=F, main=input_file)
dev.off()
