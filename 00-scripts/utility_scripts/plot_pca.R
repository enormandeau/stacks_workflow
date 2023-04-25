#!/usr/bin/env Rscript

# Usage:
#     <program> prefix popmap npc_to_plot output
#     prefix : prefix to  012 files; ex : test to read test.012, test.012.indv 7 test.012.pos
#     popmap :: 2tab delimited column: first = ind name, second = pop. no header
# 
# Note:
#     FILE SHOULD BE IMPUTED, will not be performed in this script Use 02_impute.py to perform simple imputation

# Code adapted from Florent Sylvestre, 2023-04-25

# Modules
library(data.table)
library(magrittr)
library(ggplot2)

# Parse user input
if(length(commandArgs(trailingOnly = T)) == 4) {

    ##read input file
    dat <- as.data.frame(fread(paste0(commandArgs(trailingOnly = T)[1], ".012"))[, -1])
    ind = as.data.frame(fread(paste0(commandArgs(trailingOnly = T)[1], ".012.indv"), header = F))
    popmap <- read.table(commandArgs(trailingOnly = T)[2])
    NPC_to_plot <- as.numeric(commandArgs(trailingOnly = T)[3])
    output <- commandArgs(trailingOnly = T)[4]

} else {

    cat("
Usage:
    <program> prefix popmap npc_to_plot output
    prefix : prefix to  012 files; ex : test to read test.012, test.012.indv 7 test.012.pos
    popmap :: 2tab delimited column: first = ind name, second = pop. no header

Note:
    FILE SHOULD BE IMPUTED, will not be performed in this script Use 02_impute.py to perform simple imputation

")

quit()

}

##running pca
ind <- merge(ind, popmap, sort=F)
pca <- prcomp(dat, scale=F)

###plots
PCA_plot <- function(pca, Group, pcvector=c(1, 2)) {

    sd <- pca$sdev/sum(pca$sdev) * 100

    pca$x %>%
    as.data.frame() %>%

    ggplot(aes_string(x = paste0("PC", pcvector[1]),
                      y = paste0("PC", pcvector[2]))) +

    geom_point(aes(col = Group, shape = Group)) +

    stat_ellipse(aes(col = as.factor(Group))) +

    theme_bw() +

    labs(x = paste0("PC", pcvector[1], "   (", trunc(sd[pcvector[1]] * 100)/100, "%)"),
         y = paste0("PC", pcvector[2], "   (", trunc(sd[pcvector[2]] * 100)/100, "%)"),
         col = "Group")

}

List_plot <- list()

for( i in seq(1: NPC_to_plot)[c(T, F)]) {

    if(i != NPC_to_plot) {

        List_plot[[paste(i)]] <- PCA_plot(pca, ind$V2, c(i, i+1))

    } else {

        List_plot[[paste(i)]] <- PCA_plot(pca, ind$V2, c(i, i-1))
    }

    ggsave(paste0(output, "_", i, ".png"), plot = List_plot[[paste(i)]], device = "png")
    ggsave(paste0(output, "_", i, ".pdf"), plot = List_plot[[paste(i)]], device = "pdf")

    #ggsave(paste0(output, "_", i, ".png"), plot = List_plot[[paste(i)]] +
    #   scale_color_manual(values = c("red", "blue")), device = "png")
    #ggsave(paste0(output, "_", i, ".pdf"), plot = List_plot[[paste(i)]] +
    #   scale_color_manual(values = c("red", "blue")), device = "pdf")
}
