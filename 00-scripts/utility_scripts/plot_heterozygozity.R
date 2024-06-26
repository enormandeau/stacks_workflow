#!/usr/bin/env Rscript

# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_figure = paste0(input_file, ".graph.png")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)

# Figure
png(output_figure, width=1000, height=600)
par(mar=c(12, 4, 2, 2) + 0.1)

boxplot(data[,1] ~ data[,3],
        outpch=19,
        outcol="black",
        main="Boxplot of heterozygosity per samples regrouped by populations",
        xlab="",
        xaxt="n",
        ylab="Heterozygosity by group from vcftools")

axis(1, at=1:length(unique(data[,3])), las=2, labels=unique(data[,3]))
dev.off()
