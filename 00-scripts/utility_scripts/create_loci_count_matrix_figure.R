# Clean memory
rm(list=ls())

# Global variables
proportion.used = 0.8
data.file = paste("count_matrix", proportion.used, "present.tsv", sep="_")
heatmap.file = paste("count_figure", proportion.used, "present.jpg", sep="_")
barplot.file = paste("barplot", proportion.used, "missing.jpg", sep="_")
histogram.file = paste("histogram", proportion.used, "missing.jpg", sep="_")
missing.file = paste("proportion", proportion.used, "missing.tsv", sep="_")

# Trick to speedup loading the data
temp = read.table(data.file, header=T, nrows=5)
classes = sapply(temp, class)
num.col = ncol(temp)
rm(temp)

# Load data
cat("Loading data...\n")
data = read.table(data.file, header=T, sep="\t", colClasses=classes)

# Use subset to create heatmap
cat("Creating figures...\n")
d = log(data[1:400,] + 1)

jpeg(heatmap.file, width=2000, height=2000)
    heatmap(as.matrix(d), scale="none", Rowv=NA, Colv=NA)
dev.off()

# Barplot of missing data per individual
res = c()
for (i in 1:ncol(data)) {
    res = c(res, length(data[data[,i] == 0.0,i]) / nrow(data))
}

jpeg(barplot.file, width=800, height=800)
    barplot(res, ylim=c(0, 1))
dev.off()

jpeg(histogram.file, width=800, height=800)
    hist(res, xlim=c(0, 1), breaks=seq(0, 1, 0.05))
dev.off()

# Write file with individual proportion of missing data
cat("Writing output file...\n")
missing = data.frame(indName=names(data), propMissing=res)
write.table(missing, missing.file, sep="\t")

