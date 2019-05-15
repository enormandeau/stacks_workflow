#!/usr/bin/env Rscript
# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = paste0(input_file, ".categoried")
output_image_file_1 = paste0(input_file, "_1.png")
output_image_file_2 = paste0(input_file, "_2.png")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)
d = data[,c("MedRatio", "PropHet", "PropHomRare", "Fis", "MedCovHet", "MedCovHom")]

black = "#00000011"
red = "#FF000022"
green = "#00AA0044"
blue = "#0000FF22"
purple = "#DD00AA22"

# All loci marked singleton before filters
d$Color = black

# Fis is too negative = singleton
d$Color[d$Fis < -0.1] = red
d$Color[d$Fis + d$MedRatio / 3 < 0.10] = red

# Very low Fis = diverged
d$Color[d$Fis < -0.6] = blue

# Very high Fis = low coverage and low confidence
d$Color[d$Fis - d$PropHet / 4 > 0.5] = purple

# Loci with high coverage
d$Color[d$MedCovHom > 40 | d$MedCovHet > 40] = green

# Extract bad loci infos
bad_snps = d$Color != black
all_loci = unique(gsub("_.*", "", data$ID))
bad_loci = unique(gsub("_.*", "", data$ID[bad_snps]))

# Categorize SNPs to filter loci with next script
data$Category = "singleton"
data$Category[d$Color == red] = "duplicated"
data$Category[d$Color == blue] = "diverged"
data$Category[d$Color == green] = "highcov"
data$Category[d$Color == purple] = "lowconfidence"

write.table(data[,c("Scaffold", "Position", "ID", "Category")],
            output_file, sep="\t", quote=F, row.names=F)

# Report number of duplicated loci
cat(paste0("\nLoci: ", length(bad_loci), " / ", length(all_loci), " duplicated, diverged...\n"))

# Report number of SNPs per category
report = table(data$Category)
cat("SNPs")
print(report)

# Plots
png(output_image_file_1, width=1200, height=950)
    plot(d[,1:4], pch=16, cex=1, col=d$Color)
invisible(dev.off())

png(output_image_file_2, width=1200, height=950)
    plot(d$PropHet, d$MedRatio, pch=19, cex=1.5, col=d$Color)
invisible(dev.off())
