#!/usr/bin/env Rscript

# Cleanup
rm(list=ls())

# Parse user input
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_image_file = paste0(input_file, ".png")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)

# Subset columns
d = data[,c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")]

black =  "#00000011"
red =    "#FF000022"
orange = "#FFCC0044"
yellow = "#DDDD0066"
green =  "#00AA0044"
blue =   "#0000FF22"
purple = "#DD00AA44"

# Duplicated loci
d$Color = black
d$Color[d$MedRatio < 0.3] = yellow
d$Color[d$MedRatio > 0.7] = yellow
d$Color[d$Fis + d$MedRatio / 3 - d$PropHomRare / 4 < 0.1] = orange
d$Color[d$Fis < -0.1] = red
d$Color[d$PropHet > 0.9] = blue
d$Color[d$Fis + d$PropHet / 4 > 0.7] = purple
d$Color[d$MedCovHom > 45 | d$MedCovHet > 45] = green

# Extract bad loci infos
bad_snps = d$Color != black
all_loci = unique(gsub("_.*", "", data$ID))
bad_loci = unique(gsub("_.*", "", data$ID[bad_snps]))

# TODO extract scaffold + position

cat(paste0("Found ", length(bad_loci), " duplicated loci out of ", length(all_loci), "\n"))

# Low cov (for exploration only, not removed)
d$Color[d$MedCovHom <= 8] = blue # | d$MedCovHet <= 9] = blue

# Plot
png(output_image_file, width=1200, height=950)
    plot(d[,1:4], pch=16, cex=1, col=d$Color)
dev.off()
