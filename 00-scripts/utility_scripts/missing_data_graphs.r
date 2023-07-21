# Create missing data graphs

# Get folder name
folder = read.table(".temp_graph_folder", stringsAsFactors=F)[1,1]
subfolder= "missing_data"
graph_name = paste(folder, subfolder, "missing_data.png", sep="/")

# Get data
data = read.table(paste(folder, "missing_data.tsv", sep="/"), header=T)

# Produce graph by population
png(graph_name, width=1000, height=500)

    # Increase bottom margin for population names
    par(mar=c(8,5,4,2))

    # Create basic boxplot
    boxplot(data$ProportionMissing ~ as.factor(data$Population),
            main="Proportion of missing genotypes by population",
            xlab="",
            ylab="Proportion of missing genotypes",
            ylim=c(0, 1),
            yaxp=c(0, 1, 10),
            col="lightgrey",
            outcol="transparent",
            las=2)

    # Add dotted lines at every 5%
    for (i in seq(0.05, 1, 0.1)) {
        abline(h=i, lty=3, col="lightgrey")
    }

    # Add darker lines every 10%
    for (i in seq(0, 1, 0.1)) {
        abline(h=i, lty=2, col="lightgrey")
    }

    # Add all data points
    points(as.factor(data$Population),
           data$ProportionMissing,
           pch=19,
           col="#00000022",
           cex=1.5)

dev.off()
