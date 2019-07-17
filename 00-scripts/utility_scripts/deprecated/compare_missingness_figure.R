data = read.table("missingness.tsv", header=T)
row.names(data) = names(data)
X11()
d = as.matrix(data)
for (i in nrow(data)) {
    d[i, i] = max(d)
}

png("missingness_16pops_no_empty.png", width=800, height=800)
    heatmap(as.matrix(data), scale="none", margin=c(8,8))
dev.off()

