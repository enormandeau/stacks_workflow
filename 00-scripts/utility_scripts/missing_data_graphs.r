# Create missing data graphs

# Get folder name
folder = read.table(".temp_graph_folder", stringsAsFactors=F)[1,1]
subfolder= "missing_data"
graph_name = paste(folder, subfolder, "missing_data.png", sep="/")

# Get data
data = read.table(paste(folder, "missing_data.tsv", sep="/"), header=T)

# Produce graph by population
png(graph_name, width=1000, height=500)

    # Create basic boxplot
    boxplot(data$ProportionMissing ~ data$Population,
            main="Proportion of missing genotypes by population",
            xlab="Population",
            ylab="Proportion of missing genotypes",
            ylim=c(0, 1),
            col="lightgrey",
            outcol="transparent")

    # Add all data points
    points(as.numeric(data$Population),
           data$ProportionMissing,
           pch=19,
           col="#00000022",
           cex=1.5)

dev.off()

