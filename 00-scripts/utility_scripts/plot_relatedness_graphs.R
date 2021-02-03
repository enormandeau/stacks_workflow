#!/usr/bin/env Rscript

# Parse user input
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
threshold = as.numeric(args[2])
cat(paste("Threshold:", threshold, "\n"))
output_file = paste0(input_file, ".outliers.csv")
output_figure = paste0(input_file, ".graph.png")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)

# Get data for same-same and different comparisons
same = data[data[,1] == data[,2], ]
diff = data[data[,1] != data[,2], ]
outliers = diff[diff[,3] > threshold, ]
write.table(outliers, output_file, quote=F, sep=",")

# Figure
png(output_figure, width=600, height=900)
    par(mfrow=c(3, 1))

    hist(same[,3],
         main="Samples against themselves",
         col="grey",
         breaks=seq(-1000, 1000, by=0.05),
         xlim=c(-1, 2))

    hist(diff[,3],
         main="Samples among themselves",
         col="grey",
         breaks=seq(-1000, 1000, by=0.05),
         xlim=c(-1, 2))

    hist(outliers[,3],
         main="Outlier samples",
         col="grey",
         breaks=seq(-1000, 1000, by=0.05),
         xlim=c(-1, 2))
dev.off()
