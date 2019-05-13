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

## Allele ratios are too far from 0.5
d$Color[d$MedRatio < 0.3] = yellow
d$Color[d$MedRatio > 0.7] = yellow

# Fine tune negative Fis using MedRatio and PropHomRare info
d$Color[d$Fis + d$MedRatio / 3 - d$PropHomRare / 4 < 0.1] = orange

# Fis is too negative
d$Color[d$Fis < -0.1] = red

# Diverged duplicated loci with all samples heterozygous
d$Color[d$PropHet > 0.8] = blue

# Very high Fis means low coverage and confidence on genotypes
d$Color[d$Fis + d$PropHet / 4 > 0.7] = purple

# Loci with high coverage
d$Color[d$MedCovHom > 40 | d$MedCovHet > 40] = green

# Extract bad loci infos
bad_snps = d$Color != black
all_loci = unique(gsub("_.*", "", data$ID))
bad_loci = unique(gsub("_.*", "", data$ID[bad_snps]))

# TODO extract scaffold + position

cat(paste0("Found ", length(bad_loci), " duplicated loci out of ", length(all_loci), "\n"))

# Low cov (for exploration only, not removed)
d$Color[d$MedCovHom <= 8] = blue

# Plot
png(output_image_file, width=1200, height=950)
    plot(d[,1:4], pch=16, cex=1, col=d$Color)
dev.off()
