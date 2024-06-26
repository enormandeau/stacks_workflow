#!/usr/bin/env Rscript

# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_figure = paste0(input_file, ".graph.png")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)

# Figure
if (length(unique(data[,3])) > 10) {

    png(output_figure, width=1900, height=500)

} else {

    png(output_figure, width=1000, height=400)
}

boxplot(data[,1] ~ data[,3],
        outpch=19,
        outcex=1.3,
        outcol="black",
        main="Boxplot of heterozygosity per samples by populations",
        xlab="Population",
        ylab="Heterozygosity by group from vcftools")

points(factor(data[,3]), data[,1], pch=4, col="#00000099")

dev.off()
