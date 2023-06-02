#!/usr/bin/env Rscript

# Usage:
#     <program> prefix popmap numPC output
#     prefix: prefix of 012 files (ex : test for test.012, test.012.indv 7 test.012.pos)
#           Prepared with vcftools -012
#     popmap: 2+ tab delimited columns WITH HEADER: first = Sample, second = GROUP
#           There can be as many columns as wanted after the first two
#     numPC: number of PCs to plot
#     output: Name of output file
#
# Note:
#     012 file must be created from imputed VCF

# Code adapted from Florent Sylvestre, 2023-04-25

# Modules
library(data.table)
library(magrittr)
library(ggplot2)

# Parse user input
if(length(commandArgs(trailingOnly = T)) == 4) {

    ##read input file
    dat <- as.data.frame(fread(paste0(commandArgs(trailingOnly=T)[1], ".012"))[, -1])
    ind = as.data.frame(fread(paste0(commandArgs(trailingOnly=T)[1], ".012.indv"), header=F))
    popmap <- read.table(commandArgs(trailingOnly=T)[2], header=T)
    numPC <- as.numeric(commandArgs(trailingOnly=T)[3])
    output <- commandArgs(trailingOnly=T)[4]

} else {

    cat("
Usage:
    <program> prefix popmap numPC output

    prefix: prefix of 012 files (ex : test for test.012, test.012.indv 7 test.012.pos)
        Prepared with vcftools -012
    popmap: 2+ tab delimited columns WITH HEADER: first = Sample, second = GROUP
        There can be as many columns as wanted after the first two
    numPC: number of PCs to plot
    output: Name of output file

Note:
    012 file must be created from imputed VCF
        ")

        quit()

}

##running pca
ind <- merge(ind, popmap, by.x=1, by.y="Sample", sort=F)
pca <- prcomp(dat, scale=F)

names(ind)[1] = "Sample"
write.table(cbind(ind, pca$x[, 1:20]), paste0(output, "_used_samples_with_pca.tsv"), quote=F, row.names=F, sep="\t")
#write.table(pca$x, "pc_loadings.tsv", quote=F, row.names=F, sep="\t")

###plots
PCA_plot <- function(pca, Group, pcvector=c(1, 2)) {

    sd <- pca$sdev/sum(pca$sdev) * 100

    pca$x %>%
        as.data.frame() %>%

        ggplot(aes_string(x = paste0("PC", pcvector[1]),
                          y = paste0("PC", pcvector[2]))) +

        geom_point(aes(col = Group)) +
        #geom_point(aes(col = Group, shape = Group %% 6)) +

        stat_ellipse(aes(col = as.factor(Group))) +

        theme_bw() +

        labs(x = paste0("PC", pcvector[1], "   (", trunc(sd[pcvector[1]] * 100)/100, "%)"),
             y = paste0("PC", pcvector[2], "   (", trunc(sd[pcvector[2]] * 100)/100, "%)"),
             col = "Group")

}

List_plot <- list()

for( i in seq(1: numPC)[c(T, F)]) {

    if(i != numPC) {

        List_plot[[paste(i)]] <- PCA_plot(pca, ind$Group, c(i, i+1))

    } else {

        List_plot[[paste(i)]] <- PCA_plot(pca, ind$Group, c(i, i-1))
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
