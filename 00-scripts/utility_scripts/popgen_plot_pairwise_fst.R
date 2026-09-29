#!/usr/bin/env Rscript

# Modules
library(tidyverse)

# Parse user input
args = commandArgs(trailingOnly=TRUE)
input_file = args[1] # Output from get_pairwise_fst.sh

# Read data
data = rev(read.table(input_file, header=T, sep="\t"))

# Create heatmap
ggp <- ggplot(data, aes(Pop1, Pop2)) +    # Create default ggplot2 heatmap
    # TODO Change axes names for "Population 1" and "Population 2"
    geom_tile(aes(fill = Fst)) +
    scale_y_discrete() +
    # Round numbers to 4 digits
    geom_text(aes(label = round(Fst, 4))) +
    labs(x = "Population 1", y = "Population 2")
    scale_fill_gradient(low = "#AACCEE", high = "#116699")

# Write to file
ggsave(paste0(tools::file_path_sans_ext(input_file), ".pdf"), ggp)
