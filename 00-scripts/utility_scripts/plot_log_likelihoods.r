# Importing data
data = data=read.table("06-stacks_rx/rxstacks_log_likelihoods.tsv")

# Creating figure
png("rxstacks_log_likelihoods.png", width=800, height=600)

    plot(data[,1],
         log10(data[,2]),
         type='l',
         main="Frequency of log likelihoods",
         xlab="log likelihood",
         ylab="frequency",
         xlim=c(-250, 0),
         xaxt="n"
         )

     axis(1, at=seq(-250, 0, by=10), las=2)

dev.off()

