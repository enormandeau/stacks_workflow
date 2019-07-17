# Simulate SNP datasets to evaluate the impact of Fis filtering

# Cleanup workspac
rm(list=ls())

# Global variables
num_populations = 40
num_samples = 40
num_loci = 1000
maf = 0.2

# Computing basic frequencies
q = maf
p = 1 - q
AA = p ^ 2
BB = q ^ 2
AB = 2 * p * q
expected = 2 * p * q

# Pick genotypes for one locus
all_fis = NULL
all_fis_min = NULL
all_fis_max = NULL

# Iterate over loci
for (l in seq(num_loci)) {
    fis = NULL
    fis_min = NULL
    fis_max = NULL

    # Iterate over populations
    for (p in seq(num_populations)) {
        heterozygotes = 0

        # Iterate over samples
        for (s in seq(num_samples)) {
            rand = runif(1)

            if (rand <= AB) {
                if (rand >= 0) {
                    heterozygotes = heterozygotes + 1
                }

            }
        }

        observed = heterozygotes / num_samples
        fis = c(fis, (expected - observed) / expected)
    }

    all_fis = c(all_fis, fis)

    fis_min = c(fis_min, min(as.numeric(fis)))
    all_fis_min = c(all_fis_min, fis_min)

    fis_max = c(fis_max, max(as.numeric(fis)))
    all_fis_max = c(all_fis_max, fis_max)
}

# Plot Fis distribution
par(mfrow=c(3, 1))
hist(all_fis, xlim=c(-2, 2), breaks=seq(-10, 10, 0.2), col="lightgrey")

# Plot distribution of worst Fis amont the populations for each locus
hist(all_fis_min, xlim=c(-2, 2), breaks=seq(-10, 10, 0.25), col="lightgrey")
abline(v=-0.5, lty=2)
hist(all_fis_max, xlim=c(-2, 2), breaks=seq(-10, 10, 0.25), col="lightgrey")
abline(v=0.5, lty=2)

#print(all_fis_min)
