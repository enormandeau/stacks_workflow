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
    sequence = seq(xmin, xmax, length.out = num.bins + 1)

    if (param == "presence") {
        sequence = seq(0, 1, 0.05)
        xmin = 0
        xmax = 1
    } else if (param == "numSNP") {
        sequence = seq(0, xmax + 1, 1)
        xmin = 0
        xmax = xmax + 1
    } else if (param == "allImbalance") {
        sequence = seq(10 * round((xmin - 20) / 10), 10 * round((xmax + 20) / 10), 1)
        xmin = -20
        xmax = 20
    } else if (param == "fis") {
        sequence = seq(round(xmin - 2), round(xmax + 2), 0.02)
        xmin = -1
        xmax = 1
    } else if (param == "maxDepth") {
        sequence = seq(10 * round((xmin - 20) / 10), 10 * round((xmax + 20) / 10), 10)
        xmin = 0
        xmax = 400
    } else if (param == "medDepth") {
        sequence = seq(2 * round((xmin - 10) / 2), 2 * round((xmax + 10) / 2), 2)
        xmin = 0
        xmax = 100
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
        if (pop == "global") {
            subfolder = "global"
        } else {
            subfolder = "populations"
        }

        # Produce graph
        a = data[data$Parameter == param & data$Population == pop, ]

        if (nrow(a) > 0) {
            graph_name = paste(folder, "/", subfolder , "/", param, "_", pop, ".png", sep="")

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
            dev.off()
        }
    }
}

