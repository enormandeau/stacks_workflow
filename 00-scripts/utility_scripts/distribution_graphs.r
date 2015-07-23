# Create distribution graphs

# Get folder name
folder = read.table(".temp_graph_folder", stringsAsFactors=F)[1,1]
print(folder)

# Get data
data = read.table(paste(folder, "graph_data.tsv", sep="/"), header=T)

# Iterate over parameters
for (param in levels(data$Parameter)) {
    a = data[data$Parameter == param & data$Population == "global", ]
    xmax = max(a[,3])
    xmin = min(a[,3])
    num.bins = 50
    if (param == "numSNP") { num.bins = xmax - xmin + 1 }
    #else if (param == "presence") { num.bins = 20 }
    else if (param == "allImbalance") { num.bins = 20 }
    else if (param == "medDepth") { num.bins = xmax - xmin + 1 }

    sequence = seq(xmin, xmax, length.out = num.bins + 1)

    if (param == "presence") {
        sequence = seq(0, 1, 0.05)
    } else if (param == "fis") {
        sequence = seq(xmin - 1, xmax + 1, 0.02)
        xmin = -1
        xmax = 1
    } else if (param == "genLikelihood") {
        sequence = seq(xmin - 10, xmax + 10, 5)
        xmin = -20
        xmax = 200
    } else if (param == "heterozygosity") {
        sequence = seq(0, 1, 0.01)
        xmin = 0
        xmax = 1
    } else if (param == "mafPopulation" || param == "mafGlobal") {
        sequence = seq(0, 0.5, 0.01)
        xmin = 0
        xmax = 0.5
    }

    # Iterate over populations
    for (pop in levels(data$Population)) {
        cat(param, pop, "\n") 

        # Produce graph
        a = data[data$Parameter == param & data$Population == pop, ]

        if (nrow(a) > 0) {
            graph_name = paste(folder, "/", param, "_", pop, ".png", sep="")

            # Opening graph file
            png(graph_name, width=800, height=600)

                x = a[,3]
                y = a[,4]
                hist(rep(x, y),
                     xlab = param,
                     ylab = "Count",
                     col  = "lightgrey",
                     xlim = c(xmin, xmax),
                     breaks = sequence,
                     main = paste(param, "for", pop))

                ## Do custom plots for each type of graph
                #plot(a[,3:4],
                #     pch=4,
                #     col="#00000088",
                #     cex=0.6,
                #     main = paste(param, "for", pop))
                #lines(lowess(a[,4] ~ a[,3], f=0.2), lwd=2, lty=1, col="red")

            dev.off()
        }
    }
}

