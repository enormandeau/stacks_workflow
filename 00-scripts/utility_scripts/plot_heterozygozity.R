#!/usr/bin/env Rscript

# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_figure = paste0(input_file, ".graph.png")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)

# Figure
png(output_figure, width=1000, height=400)

boxplot(data[,1] ~ data[,3],
        outpch=19,
        outcol="black",
        main="Boxplot of heterozygosity per samples regrouped by populations",
        xlab="Population",
        ylab="Heterozygosity from vcftools")
dev.off()
