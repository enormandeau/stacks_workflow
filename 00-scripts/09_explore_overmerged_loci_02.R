#!/usr/bin/env Rscript

# Cleanup
rm(list=ls())

# Load data
data = read.table("filtered_m4_p70.overmerged.infos", header=T, stringsAsFactors=F)

# Subset columns
d = data[,c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")]

black =  "#00000011"
red =    "#FF000022"
orange = "#FFCC0066"
yellow = "#DDFF0022"
green =  "#00FF6688"
blue =   "#0000FF22"
purple = "#FF00AA22"

# Duplicated loci
d$Color = black
d$Color[d$MedRatio < 0.3] = yellow
d$Color[d$MedRatio > 0.7] = yellow
d$Color[d$PropHet > 0.55] = orange
d$Color[d$Fis + d$MedRatio / 4 < 0.05] = orange
d$Color[d$Fis < -0.1] = red
d$Color[d$MedCovHom > 50 | d$MedCovHet > 50] = green

# Extract bad loci infos
bad_snps = d$Color != black
all_loci = unique(gsub("_.*", "", data$ID))
bad_loci = unique(gsub("_.*", "", data$ID[bad_snps]))

# TODO extract scaffold + position

cat(paste0("Found ", length(bad_loci), " duplicated loci out of ", length(all_loci), "\n"))

# Low cov and too many heterozygotes
d$Color[d$Fis + d$PropHet / 4 > 0.4] = purple
d$Color[d$MedCovHom <= 10] = blue # | d$MedCovHet <= 9] = blue

# Plot
png("duplicated_loci.png", width=1200, height=950)
    plot(d[,1:4], pch=16, cex=1, col=d$Color)
dev.off()
