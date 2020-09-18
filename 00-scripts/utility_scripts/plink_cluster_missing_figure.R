#!/usr/bin/env Rscript

# Parse user input
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file1 = paste0(input_file, "_missing_1.png")
output_file2 = paste0(input_file, "_missing_2.png")

# Load data
data = read.table(input_file, header=T, stringsAsFactors=F)

print(head(data))

# Figure 1
png(output_file1, width=2000, height=2000)

minx = min(data$C1)
maxx = max(data$C1)
miny = min(data$C2)
maxy = max(data$C2)

plot(1, 1,
     main="TITLE",
     xlab="C1",
     ylab="C2",
     xlim=c(minx, maxx),
     ylim=c(miny, maxy),
     type='n')

text(data$C1, data$C2, labels=data$FID)
dev.off()

# Figure 2
png(output_file2, width=2000, height=2000)

minx = min(data$C3)
maxx = max(data$C3)
miny = min(data$C4)
maxy = max(data$C4)

plot(1, 1,
     main="TITLE",
     xlab="C3",
     ylab="C4",
     xlim=c(minx, maxx),
     ylim=c(miny, maxy),
     type='n')

text(data$C3, data$C4, labels=data$FID)
dev.off()
