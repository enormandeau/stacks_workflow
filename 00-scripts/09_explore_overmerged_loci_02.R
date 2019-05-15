#!/usr/bin/env Rscript
# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = paste0(input_file, ".categorized")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)
d = data[,c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")]

singleton =     "#00000011" # black
duplicated =    "#FF000022" # red
diverged =      "#0000FF22" # blue
lowconfidence = "#DD00AA22" # purple
highcov =       "#00AA0022" # green
#test =          "#FFAA0022" # orange

# All loci marked singleton before filters
d$Color = singleton

# MedRatio is high/low and at least one rare allele homozygote
d$Color[d$MedRatio < 0.40 & d$PropHomRare > 0.00] = lowconfidence
d$Color[d$MedRatio > 0.60 & d$PropHomRare > 0.00] = lowconfidence

# Fis is too negative = duplicated
d$Color[d$Fis < -0.1] = duplicated
d$Color[d$Fis + d$MedRatio / 3 < 0.11] = duplicated

# Very low Fis = diverged
d$Color[d$Fis < -0.6] = diverged

# Very high Fis = low coverage and low confidence
#d$Color[d$Fis - d$PropHet / 4 > 0.45] = lowconfidence
d$Color[d$Fis > 0.49] = lowconfidence

# Loci with high coverage
d$Color[d$MedCovHom > 40 | d$MedCovHet > 40] = highcov

# Extract bad loci infos
bad_snps = d$Color != singleton
all_loci = unique(gsub("_.*", "", data$ID))
bad_loci = unique(gsub("_.*", "", data$ID[bad_snps]))

# Categorize SNPs to filter loci with next script
data$Category = "singleton"
data$Category[d$Color == duplicated] = "duplicated"
#data$Category[d$Color == test] = "duplicated"
data$Category[d$Color == diverged] = "diverged"
data$Category[d$Color == lowconfidence] = "lowconfidence"
data$Category[d$Color == highcov] = "highcov"

write.table(data[,c("Scaffold", "Position", "ID", "Category")],
            output_file, sep="\t", quote=F, row.names=F)

# Report number of duplicated loci
cat(paste0("\nLoci: ", length(bad_loci), " / ", length(all_loci), " duplicated, diverged...\n"))

# Report number of SNPs per category
report = table(data$Category)
cat("SNPs")
print(report)

# Plots
png(paste0(input_file, "_1.png"), width=1200, height=950)
    plot(d[,1:4], pch=16, cex=1, col=d$Color)
invisible(dev.off())

png(paste0(input_file, "_2.png"), width=1200, height=950)
    plot(d$PropHet, d$MedRatio, pch=19, cex=1.5, col=d$Color, xlim=c(0, 1), ylim=c(0, 0.8))
invisible(dev.off())

single = d[data$Category == "singleton", ]# | data$Category == "highcov" | data$Category == "test", ]
png(paste0(input_file, "_3.png"), width=1200, height=950)
    plot(single$PropHet,
         single$MedRatio,
         pch=19, cex=1.5, col=single$Color, xlim=c(0, 1), ylim=c(0, 0.8))
invisible(dev.off())
