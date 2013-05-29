#!/usr/bin/python
"""Filter sumstats.tsv files from STACKS to remove :

Usage:
  ./filterStacksSNPs.py inFile maxSnpNumber maxAlleleNumber minPresence maxHetero minAlleleFreq minFis maxFis

inFile = batch_1.sumstat.tsv or similarly names output of STACKS (v0.99995+)
maxSnpNumber = maximum number of SNPs in a single locus (int, 1 or more)
maxAlleleNumber = maximum number of possible alleles per SNP (int, 2 to 4)
minPresence = minimal number of present individuals (int, 1 or more)
maxHetero = maximal frequency of heterozygous individuals (float, 0 to 1)
minAlleleFreq = minimum frequency of rare allele (float, 0 to 1)
minFis = minimum allowed Fis value (float, -1 to 1)
maxFis = maximum allowed Fis value (float, -1 to 1)
"""

# Importing modules
import sys
from collections import defaultdict

# Defining classes
class SNP(object):
    """Store and manipulate locus information
    """
    loci = defaultdict(lambda: defaultdict(dict))
    def __init__(self, line):
        l = line.strip().split("\t")
        self.line = l
        self.locus = l[1]
        self.position = int(l[4])
        self.population = l[5]
        self.maxFreq = float(l[9])
        self.minFreq = 1 - self.maxFreq
        self.obsHet = float(l[10])
        self.expHet = float(l[12]) 
        self.fis = float(l[17])
        try:
            self.presence = int(l[8])
        except:
            self.presence = 0.0
        try:
            self.hetero = float(l[10])
        except:
            self.hetero = 0.0
        SNP.loci[self.locus][self.position][self.population] = self
    def __str__(self):
        return "\t".join(self.line) + "\n"
            
# Conditional functions
def filter_true(lst): # filter true results
    return [l for l in lst if l]

def at_least_n_true(lst, n):
    return len(filter_true(lst)) >= n

# Filter functions
def max_snp_number(dct, n):
    for locus, positions in dct.items():
        if len(positions) > n:
            dct.pop(locus)

def enough_individuals(snp, n):
    return snp.presence >= n

def max_heterozygote_freq(snp, freq):
    return snp.hetero <= freq

def min_allele_freq(snp, freq):
    return snp.minFreq >= freq

def min_fis(snp, lower):
    return snp.fis >= lower

def max_fis(snp, upper):
    return snp.fis <= upper

# Trim workhouse
def trim_snp_dict(dct, fun, arg, cond):
    """dct = SNP.loci
       fun = test function (eg: enough_individuals)
       arg = argument to fun. If none needed, put None
       cond = conditional function (eg: at_least_some_true)
    """
    for locus, positions in dct.items():
        for position, snp in positions.items():
            if not cond([fun(x, arg) for x in snp.values()]):
                if locus in dct:
                    dct.pop(locus)

# Remove SNPs that are not present in both populations
# TODO Generalize for more than 2 populations
def all_pops(dct):
    for locus, positions in dct.items():
        for position, population in positions.items():
            if len(dct[locus][position]) < 2:
                dct[locus].pop(position)

# Add an Fst to the line of info of each SNP
# TODO Generalize for more than 2 populations
def add_fst(dct):
    #print "n1 n2 n_tot p1 p2 p_tot He_tot He_mean Fst"
    for locus, positions in dct.items():
        for position, population in positions.items():
            n1 = population["1"].presence
            n2 = population["2"].presence
            n_tot = n1 + n2
            p1 = population["1"].maxFreq
            p2 = population["2"].maxFreq
            p_tot = (p1 * n1 + p2 * n2) / n_tot # Useful to normalize
            He_tot = 2 * p_tot * (1 - p_tot)
            He1 = population["1"].expHet
            He2 = population["2"].expHet
            He_mean = (He1 * n1 + He2 * n2) / n_tot
            try:
                Fst = (He_tot - He_mean) / He_tot
            except:
                Fst = 0
            population["1"].line.append(str(He_tot))
            population["1"].line.append(str(Fst))
            population["2"].line.append(str(He_tot))
            population["2"].line.append(str(Fst))
            #print n1, n2, n_tot, p1, p2, p_tot, He_tot, He_mean, Fst
            
# Counting remaining snps
def count_snps(dct):
    return sum([len(x) for x in dct.values()])

# Write good SNPs to file
def write_remaining_snps(filename, snp_dict):
    with open(filename, "w") as outf:
        for h in output_header:
            outf.write(h)
        for locus, positions in snp_dict:
            for position, populations in positions.items():
                for population, snp in populations.items():
                    outf.write(str(snp))

# Main function
if __name__ == "__main__":
    try: # Parse user input
        inputFile = sys.argv[1]
        maxSnpNumber = int(sys.argv[2])
        maxAlleleNumber = int(sys.argv[3])
        minPresence = int(sys.argv[4])
        maxHetero = float(sys.argv[5])
        minAlleleFreq = float(sys.argv[6])
        minFis = float(sys.argv[7])
        maxFis = float(sys.argv[8])
    except:
        print __doc__
        sys.exit(1)
    
    # Asserting that parmeters have sensible values
    try:
        with open(inputFile) as test:
            pass
    except:
        print "Input error: file not found (" + inputFile + ")"
        sys.exit(1)

    assert maxSnpNumber > 0, "maxSnpNumber must be a non-null integer"
    assert maxAlleleNumber in [2, 3, 4], "maxAlleleNumber must be 2, 3 or 4"
    assert minPresence >= 1, "minPresence must be a a non-number integer"
    assert maxHetero > 0 and maxHetero <= 1, "maxHetero must be a decimal between 0 and 1"
    assert minAlleleFreq > 0 and minAlleleFreq <= 0.5, "minAlleleFreq must be a decimal between 0 and 0.5"
    assert minFis >= -1 and minFis <= 1, "minFis must be a decimal between -1 and 1"
    assert maxFis >= -1 and maxFis <= 1, "maxFis must be a decimal between -1 and 1"

    # Constructing the SNP instances and dictionary
    output_header = ""
    with open(inputFile) as f:
        for line in f:
            if not line.startswith("#"):
                SNP(line) 
            else:
                output_header = line.strip() + "\tHe_tot\tFst\n"

    # Initialize blacklist
    blacklist = set()
    for s in SNP.loci.keys():
        blacklist.add(s)

    # Creating output file header
    header = ["File",
              "1orMoreSNP",
              "AtMost{0}SNP".format(maxSnpNumber),
              "Biallelic",
              "AtLeast{0}Ind".format(minPresence),
              "Max{0}Hetero".format(maxHetero),
              "Min{0}AlleleFreq".format(minAlleleFreq),
              "minFis{}".format(minFis),
              "maxFis{}".format(maxFis)]

    # Initializing the filtering output
    to_print = [inputFile]

    ### Start filtering 
    # Remove SNPs that are not present in both populations
    all_pops(SNP.loci)

    # Add calculated Fst value at the end of the line
    add_fst(SNP.loci)

    # >= 1 SNP: no filter
    to_print.append(str(count_snps(SNP.loci)))
    write_remaining_snps("filtered_loci_1_allSNPs.tsv", SNP.loci.items())

    # <= maxSnpNumber SNP
    max_snp_number(SNP.loci, maxSnpNumber)
    to_print.append(str(count_snps(SNP.loci)))
    write_remaining_snps("filtered_loci_2_maxnSNPs.tsv", SNP.loci.items())
    
    # <= maxAlleleNumber alleles per SNP
    # TODO This filter not implemented
    to_print.append(str(count_snps(SNP.loci)))
    write_remaining_snps("filtered_loci_3_maxAlleleNum.tsv", SNP.loci.items())

    # Loop over a list of filters [test, parameter, conditional, testName]
    trim_step = 4
    trim_jobs = [[enough_individuals, minPresence, all, header[4]], 
                 [max_heterozygote_freq, maxHetero, all, header[5]],
                 [min_allele_freq, minAlleleFreq, any, header[6]],
                 [min_fis, minFis, all, header[7]],
                 [max_fis, maxFis, all, header[8]]]
    
    for job in trim_jobs:
        trim_snp_dict(SNP.loci, job[0], job[1], job[2])
        to_print.append(str(count_snps(SNP.loci)))
        (write_remaining_snps("filtered_loci_" + str(trim_step) +
            "_" + job[3] + ".tsv", SNP.loci.items()))
        trim_step += 1

    # Update whitelist and blacklist 
    whitelist = set()
    for s in SNP.loci.keys():
        whitelist.add(s)
    blacklist = blacklist.difference(whitelist)

    # Print whitelist and blacklist to files
    with open("whitelist.txt", "w") as f:
        for l in sorted(list(whitelist)):
            f.write(str(l) + "\n")

    with open("blacklist.txt", "w") as f:
        for l in sorted(list(blacklist)):
            f.write(str(l) + "\n")

    # Print results to screen
    print "\t".join(header)
    print "\t".join(to_print)

