#!/usr/bin/env Rscript
# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = paste0(input_file, ".categorized")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)
d = data[,c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHom", "MedCovHet")]
for (column in names(d)) {
    d[[column]][!is.finite(d[[column]])] = NA_real_
}

canonical =     "#00000011" # black
duplicated =    "#FF000022" # red
diverged =      "#0000FF22" # blue
lowconf =       "#DD00BB22" # purple
highcov =       "#00AA0022" # green
mas =           "#FFAA0022" # orange
insufficient =  "#77777788" # grey

# All loci marked canonical before filters
d$Color = canonical
# Canonical requires usable ratio and Fis evidence. Independent criteria below
# may still assign a category when their own required statistics are available.
d$Color[which(!is.finite(d$MedRatio) | !is.finite(d$Fis))] = insufficient

# Loci with high coverage
maxMedCovHom = 50
maxMedCovHet = 200
d$Color[which(d$MedCovHom > maxMedCovHom | d$MedCovHet > maxMedCovHet)] = highcov
d$MedCovHom[which(d$MedCovHom > maxMedCovHom * 2)] = maxMedCovHom * 2
d$MedCovHet[which(d$MedCovHet > maxMedCovHet * 2)] = maxMedCovHet * 2

# MedRatio is high/low and at least one rare allele homozygote
d$Color[which(d$MedRatio < 0.20)] = lowconf # & d$PropHomRare > 0.00] = lowconf
d$Color[which(d$MedRatio > 0.80)] = lowconf # & d$PropHomRare > 0.00] = lowconf

# Fis is too negative = duplicated
d$Color[which(d$Fis < -0.4)] = duplicated
d$Color[which(d$Fis + d$MedRatio < 0.08)] = duplicated
d$Color[which(d$Fis + d$MedRatio * 3 < 0.78)] = duplicated
d$Color[which(d$Fis + d$MedRatio * 8 < 2.3)] = duplicated

# Very low Fis = diverged
d$Color[which(d$Fis < -0.8)] = diverged
d$Color[which(d$Fis + d$MedRatio * 2 < -0.00)] = diverged
d$Color[which(d$Fis + d$MedRatio * 3 < 0.20)] = diverged
d$Color[which(d$Fis + d$MedRatio * 8 < 1.5)] = diverged

# High Fis
d$Color[which(d$Fis > 0.9)] = lowconf

# Too few samples with rare allele
# Do not interpret an entirely uncalled site as evidence for a rare allele.
d$Color[which(is.finite(d$PropHet) & data$NumHet + data$NumRare < 3)] = mas

# Extract bad loci infos
bad_snps = d$Color != canonical
all_loci = unique(gsub("_.*", "", data$ID))
bad_loci = unique(gsub("_.*", "", data$ID[bad_snps]))

# Categorize SNPs to filter loci with next script
data$Category = "canonical"
data$Category[d$Color == duplicated] = "duplicated"
data$Category[d$Color == mas] = "mas"
data$Category[d$Color == diverged] = "diverged"
data$Category[d$Color == lowconf] = "lowconf"
data$Category[d$Color == highcov] = "highcov"
data$Category[d$Color == insufficient] = "insufficient_data"

write.table(data[,c("Scaffold", "Position", "ID", "Category")],
            output_file, sep="\t", quote=F, row.names=F)

# Report number of SNPs per category
report = table(data$Category)
cat("SNPs")
print(report)

# Plots
png(paste0(input_file, "_1.png"), width=1600, height=1150)
    available = vapply(d[,1:6], function(x) any(is.finite(x)), logical(1))
    if (sum(available) >= 2 && nrow(d) > 0) {
        plot(d[,which(available),drop=FALSE], pch=16, cex=0.6, col=d$Color)
    } else {
        plot.new()
        text(0.5, 0.5, "Insufficient finite statistics for pairwise plots")
    }
invisible(dev.off())

png(paste0(input_file, "_2.png"), width=1600, height=1150)
    plot(d$PropHet, d$MedRatio, pch=19, cex=1.5, col=d$Color, xlim=c(0, 1), ylim=c(0, 0.8))
invisible(dev.off())

single = d[data$Category == "canonical", ]
png(paste0(input_file, "_3.png"), width=1600, height=1150)
    plot(single$PropHet,
         single$MedRatio,
         pch=19, cex=1.5, col=single$Color, xlim=c(0, 1), ylim=c(0, 0.8))
invisible(dev.off())
