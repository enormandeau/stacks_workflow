#!/usr/bin/env Rscript
# Code adapted from Florent Sylvestre, 2023-04-25

# Modules
library(data.table)
library(magrittr)
library(ggplot2)

# Parse user input
if(length(commandArgs(trailingOnly = T)) == 6) {

    ##read input file
    dat <- as.data.frame(fread(paste0(commandArgs(trailingOnly=T)[1], ".012"))[, -1])
    ind = as.data.frame(fread(paste0(commandArgs(trailingOnly=T)[1], ".012.indv"), header=F))
    popmap <- read.table(commandArgs(trailingOnly=T)[2], header=T)
    numPC <- as.numeric(commandArgs(trailingOnly=T)[3])
    variable1 <- as.character(commandArgs(trailingOnly=T)[4])
    variable2 <- as.character(commandArgs(trailingOnly=T)[5])
    output <- commandArgs(trailingOnly=T)[6]

    cat(paste0("Producing figure for: ", variable1, " and ", variable2, "\n"))

} else {

    cat("
Usage:
    <program> prefix popmap numPC variable1 variable2 output

    prefix: prefix of 012 files (ex : test for test.012, test.012.indv 7 test.012.pos)
          Prepared with vcftools -012
    popmap: 2+ tab delimited columns WITH HEADER: first = Sample
          There can be as many columns as wanted after the first two but two of them
          must be the `variable1` and `variable2` of interest
    numPC: number of PCs to plot
    variable1: Name of variable from popmap file
    variable2: Name of variable from popmap file
    output: Name of output file

Note:
    012 file must be created from imputed VCF
        ")

        quit()

}

##running pca
ind <- merge(ind, popmap, by.x=1, by.y="Sample", sort=F)
ind[, variable1] = as.factor(ind[, variable1])
ind[, variable2] = as.factor(ind[, variable2])

pca <- prcomp(dat, scale=F)

names(ind)[1] = "Sample"
write.table(cbind(ind, pca$x[, 1:20]), paste0(output, "_used_samples_with_pca.tsv"), quote=F, row.names=F, sep="\t")
write.table(pca$x, "pc_loadings.tsv", quote=F, row.names=F, sep="\t")

###plots
PCA_plot <- function(pca, ind, v1, v2, pcvector=c(1, 2)) {

    sd <- pca$sdev/sum(pca$sdev) * 100

    pca$x %>%
        as.data.frame() %>%

        ggplot(aes_string(x = paste0("PC", pcvector[1]),
                          y = paste0("PC", pcvector[2]))) +

        geom_point(aes(col = ind[, v1], shape = ind[, v2]), size=3, alpha=0.6, fill=NA) +

        theme_bw() +

        labs(x = paste0("PC", pcvector[1], "   (", trunc(sd[pcvector[1]] * 100)/100, "%)"),
             y = paste0("PC", pcvector[2], "   (", trunc(sd[pcvector[2]] * 100)/100, "%)"),
             col = v1,
             shape = v2)
}

List_plot <- list()

for( i in seq(1: numPC)[c(T, F)]) {

    if(i != numPC) {

        List_plot[[paste(i)]] <- PCA_plot(pca, ind, variable1, variable2, c(i, i+1))

    } else {

        List_plot[[paste(i)]] <- PCA_plot(pca, ind, variable1, variable2, c(i, i+1))
    }

    ggsave(paste0(output, "_", (i+1)/2, ".png"),
           plot = List_plot[[paste(i)]],
           device = "png",
           height=8,
           width=10,
           units="in")

    ggsave(paste0(output, "_", (i+1)/2, ".pdf"),
           plot = List_plot[[paste(i)]],
           device = "pdf",
           height=8,
           width=10,
           units="in")
}
