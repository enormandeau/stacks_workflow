#!/usr/bin/env Rscript
# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = paste0(input_file, ".categorized")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)
d = data[,c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")]

singleton =     "#00000005" # black
duplicated =    "#FF000022" # red
diverged =      "#0000FF22" # blue
lowconf =       "#DD00AA22" # purple
highcov =       "#00AA0022" # green
mas =           "#FFAA0022" # orange

# All loci marked singleton before filters
d$Color = singleton

# MedRatio is high/low and at least one rare allele homozygote
d$Color[d$MedRatio < 0.3] = lowconf # & d$PropHomRare > 0.00] = lowconf
d$Color[d$MedRatio > 0.7] = lowconf # & d$PropHomRare > 0.00] = lowconf

# Fis is too negative = duplicated
d$Color[d$Fis < -0.1] = duplicated
#d$Color[d$Fis + d$MedRatio < 0.08] = duplicated
#d$Color[d$Fis + d$MedRatio * 3 < 0.78] = duplicated
#d$Color[d$Fis + d$MedRatio * 8 < 2.3] = duplicated

# Very low Fis = diverged
d$Color[d$Fis < -0.8] = diverged
#d$Color[d$Fis + d$MedRatio * 2 < -0.00] = diverged
#d$Color[d$Fis + d$MedRatio * 3 < 0.20] = diverged
#d$Color[d$Fis + d$MedRatio * 8 < 1.5] = diverged

# High Fis
d$Color[d$Fis > 0.90] = lowconf

# Loci with high coverage
d$Color[d$MedCovHom > 60 | d$MedCovHet > 60] = highcov

# Too few samples with rare allele
d$Color[data$NumHet + data$NumRare < 3] = mas

# Extract bad loci infos
bad_snps = d$Color != singleton
all_loci = unique(gsub("_.*", "", data$ID))
bad_loci = unique(gsub("_.*", "", data$ID[bad_snps]))

# Categorize SNPs to filter loci with next script
data$Category = "singleton"
data$Category[d$Color == duplicated] = "duplicated"
data$Category[d$Color == mas] = "mas"
data$Category[d$Color == diverged] = "diverged"
data$Category[d$Color == lowconf] = "lowconf"
data$Category[d$Color == highcov] = "highcov"

write.table(data[,c("Scaffold", "Position", "ID", "Category")],
            output_file, sep="\t", quote=F, row.names=F)

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

single = d[data$Category == "singleton", ]
png(paste0(input_file, "_3.png"), width=1200, height=950)
    plot(single$PropHet,
         single$MedRatio,
         pch=19, cex=1.5, col=single$Color, xlim=c(0, 1), ylim=c(0, 0.8))
invisible(dev.off())
