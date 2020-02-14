#!/usr/bin/env Rscript
# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = paste0(input_file, ".categorized")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)
#d = data[,c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")]
d = data[,c("AvgRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")]
names(d) = c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")

singleton =     "#00000011" # black
duplicated =    "#FF000022" # red
diverged =      "#0000FF22" # blue
lowconf =       "#DD00DD44" # purple
highcov =       "#00AA0088" # green
mas =           "#FFAA0044" # orange

# All loci marked singleton before filters
d$Color = singleton

# Fis is too negative = duplicated
d$Color[d$Fis < -0.30] = duplicated
d$Color[d$MedRatio > 0.22 & d$Fis < -0.45] = duplicated
d$Color[d$MedRatio > 0.22 & d$Fis + d$MedRatio < 0.10] = duplicated
d$Color[d$MedRatio > 0.22 & d$Fis + d$MedRatio * 3 < 1.00] = duplicated
d$Color[d$MedRatio > 0.22 & d$Fis + d$MedRatio * 8 < 2.85] = duplicated

# Very low Fis = diverged
d$Color[d$Fis < -0.1 & d$Fis + d$MedRatio < -0.175] = diverged
d$Color[d$Fis < -0.1 & d$Fis + d$MedRatio * 3 < 0.54] = diverged
d$Color[d$Fis < -0.1 & d$Fis + d$MedRatio * 8 < 1.70] = diverged
d$Color[d$Fis < -0.1 & d$MedRatio < 0.22] = diverged

# Salvage back some singletons
d$Color[d$Fis > 0.0 & d$MedRatio + d$Fis * 2 > 0.35 & d$MedRatio < 0.40] = singleton

# MedRatio is high/low and Fis is too positive
d$Color[d$Fis > 0.0 & d$MedRatio + d$Fis * 2 > 0.35 & d$MedRatio < 0.28] = lowconf
d$Color[d$Fis > 0.0 & d$MedRatio + d$Fis * 2 > 0.35 & d$MedRatio > 0.7] = lowconf
d$Color[d$Fis > -0.1 & d$MedRatio < 0.22] = lowconf
d$Color[d$Fis > -0.1 & d$MedRatio > 0.70] = lowconf

# Loci with high coverage
d$Color[d$MedCovHom > 100 | d$MedCovHet > 100] = highcov

# Very low Fis = diverged
d$Color[d$Fis < -0.65] = diverged

# Too few samples with rare allele
d$Color[data$NumHet + data$NumRare < 8] = mas

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
cat("SNPs")
print(table(data$Category))

# Plots
#png(paste0(input_file, "_1.png"), width=1200, height=950)
png(paste0(input_file, "_1.png"), width=1200, height=950)
    plot(d[,1:4], pch=16, cex=0.8, col=d$Color)
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
