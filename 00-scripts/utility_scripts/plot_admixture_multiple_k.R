#!/usr/bin/Rscript

# TODO Permit list of specific K values separated by commas

library(RColorBrewer)
color_scheme="Spectral"

# Usage: plotADMIXTURE.r -p <prefix>
#                       -i <info file, 2-column file with header and \t pop>
#                       -k <max K value>
#                       -l <comma-separated list of populations/species in the order to be plotted>

# This R script makes barplots for K=2 and all other K values until max K
# (specified with -k). It labels the individuals and splits them into
# populations or species according to the individual and population/species
# names in the 2-column file specified with -i.  The order of
# populations/species follows the list of populations/species given with -l.

# Usage example: plotADMIXTURE.r -p fileXY -i file.ind.pop.txt -k 4 -pop pop1,pop2,pop3

# In this example, the script would use the files fileXY.2.Q, fileXY.3.Q, fileXY.4.Q to make barplots for the three populations.
# file.ind.pop.txt should contain one line for each individual in the same order as in the admixture files e.g.
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop3

# Author: Joana Meier, September 2019
# Modified by Eric Normandeau, 2023-04-24

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="prefix name (with path if not in the current directory)", metavar="character"),
  make_option(c("-i", "--infofile"), type="character", default=NULL,
              help="info text file containing for each individual the population/species information", metavar="character"),
  make_option(c("-k", "--maxK"), type="integer", default=NULL,
              help="maximum K value", metavar="integer"),
  make_option(c("-m", "--minK"), type="integer", default=2,
              help="minimum K value", metavar="integer"),
  make_option(c("-l", "--populations"), type="character", default=NULL,
              help="comma-separated list of populations/species in the order to be plotted", metavar="character"),
  make_option(c("-o", "--outPrefix"), type="character", default="default",
              help="output prefix (default: name provided with prefix)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide the prefix", call.=FALSE)
} else if (is.null(opt$infofile)){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
} else if (is.null(opt$maxK)){
  print_help(opt_parser)
  stop("Please provide the maximum K value to plot", call.=FALSE)
} else if (is.null(opt$populations)){
  print_help(opt_parser)
  stop("Please provide a comma-separated list of populations/species", call.=FALSE)
}

# If no output prefix is given, use the input prefix
if(opt$outPrefix=="default") opt$outPrefix=opt$prefix

# Assign the first argument to prefix
prefix=opt$prefix

# Get individual names in the correct order
labels = read.table(opt$infofile, header=T)[, c(1, 2)]

# Name the columns
names(labels) = c("ind", "pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n = factor(labels$pop, levels=unlist(strsplit(opt$populations, ",")))
levels(labels$n) = c(1: length(levels(labels$n)))
labels$n = as.integer(as.character(labels$n))

# read in the different admixture output files
minK=opt$minK
maxK=opt$maxK
tbl = lapply(minK:maxK, function(x) read.table(paste0(prefix, ".", x, ".Q")))

# Prepare spaces to separate the populations/species
rep = as.vector(table(labels$n))
spaces = 0
inter_pop_space = 1

for(i in 1:length(rep)) {
    spaces=c(spaces, rep(0, rep[i]-1), inter_pop_space)
}

spaces = spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual for each K as a separate row
#png(file=paste0(opt$outPrefix, ".png"), width = 3000, height = 1800, res=200)
pdf(file=paste0(opt$outPrefix, ".pdf"), width = 8, height = 5)

    par(mfrow=c(maxK - 1, 1),
        mar=c(0, 1, 0, 0),
        oma=c(2, 1, 1, 1),
        mgp=c(0, 0.2, 0),
        xaxs="i",
        cex.lab=1.2,
        cex.axis=0.8)

    # Plot minK
    bp = barplot(t(as.matrix(tbl[[1]][order(labels$n), ])),
              col=brewer.pal(n=minK, "Spectral"),
              xaxt="n",
              border=NA,
              ylab=paste0("K=", minK),
              yaxt="n",
              space=spaces)

    # Add sample names at top
    # axis(3, at=bp, labels=labels$ind[order(labels$n)], las=2, tick=F, cex=0.6, )
    # Add to par(mfrow above to get sample names
    # oma=c(2, 1, 9, 1),

    # Plot higher K values
    if(maxK > minK) {

        lapply(2: (maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n), ])),
                                                col=brewer.pal(n=x+1, "Spectral"),
                                                xaxt="n",
                                                border=NA,
                                                ylab=paste0("K=", x+1),
                                                yaxt="n",
                                                space=spaces))
    }

    axis(1, at=c(which(spaces==inter_pop_space),
         bp[length(bp)])-diff(c(1, which(spaces==inter_pop_space),bp[length(bp)]))/2,
         labels=unlist(strsplit(opt$populations, ",")), cex=2)

dev.off()
