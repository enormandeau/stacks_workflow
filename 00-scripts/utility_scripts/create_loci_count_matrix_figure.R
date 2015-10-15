# Clean memory
rm(list=ls())

# Global variables
data.file = "count_matrix.tsv"

# Trick to speedup loading the data
temp = read.table(data.file, header=T, nrows=5)
classes = sapply(temp, class)
num.col = ncol(temp)
rm(temp)

# Load data
data = read.table(data.file, header=T, sep="\t", colClasses=classes)

# Use subset to create heatmap
d = data[1:400, ]

jpeg("count_figure_all_present.jpg")
    heatmap(as.matrix(d))
dev.off()

