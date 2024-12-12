#!/usr/bin/env python3
"""From files created by extracting the first 8 columns of non-commented lines
in two VCFs and a reference, format sequences for panel design purposes.

Usage:
    <program> input_genome input_complete_info input_filtered_info min_maf flank_size output_sequences
"""

# Modules
from collections import defaultdict
from random import choice
import gzip
import sys

# Classes
class Fasta(object):
    """Fasta object with name and sequence
    """

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

    def __repr__(self):
        return self.name + " " + self.sequence[:31]

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator
    """
    with myopen(input_file) as f:
        sequence = []
        name = ""
        begun = False

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if begun:
                    yield Fasta(name, "".join(sequence))

                name = line[1:]
                sequence = ""
                begun = True

            else:
                sequence += line

        if name != "":
            yield Fasta(name, "".join(sequence))

def compression_ratio(seq, minimum_ratio=0.0, maximum_ratio=1.0):
    ratio = len(gzip.compress(seq.upper().encode())) / len(seq)
    return (ratio - minimum_ratio) / (maximum_ratio - minimum_ratio)

# Parse user input
try:
    input_genome = sys.argv[1]
    input_complete_info = sys.argv[2]
    input_filtered_info = sys.argv[3]
    min_maf = float(sys.argv[4])
    flank_size = int(sys.argv[5])
    output_sequences = sys.argv[6]
except:
    print(__doc__)
    sys.exit(1)

# Load all SNPs in dict by regions
all_snps = defaultdict(lambda: defaultdict(list))
kept_snps = 0

with open(input_complete_info, "rt") as infile:
    for line in infile:
        chrom, pos, locus, allele1, allele2, _, _, infos = line.strip().split()
        pos = int(pos)
        locus_num, locus_pos, locus_sens = locus.split(":")
        num_samples, maf = infos.split(";")
        num_samples = int(num_samples.split("=")[1])
        maf = float(maf.split("=")[1])

        snp = (chrom, pos, locus_num, locus_pos, locus_sens, allele1, allele2, num_samples, maf)

        if maf >= min_maf:
            kept_snps += 1
            all_snps[chrom][pos // flank_size].append(snp)

print(f"Kept {kept_snps} SNPs above MAF={min_maf}")

# Load genome
chromosomes = dict()
sequences = fasta_iterator(input_genome)

for s in sequences:
    chromosomes[s.name] = s.sequence

print(f"Loaded {len(chromosomes)} sequences from {input_genome}")

# Normalize complexity values by maximum for window size
# Generate random sequences of flank_size and find maximum_ratio
minimum_ratio = compression_ratio("".join("A" * (2*flank_size+1)))
maximum_ratio = 0.0

for i in range(50):
    random_sequence = "".join(choice("ACGT") for _ in range(2*flank_size+1))
    ratio = compression_ratio(random_sequence)
    maximum_ratio = max(ratio, maximum_ratio)

# Read filtered SNPs and write sequences for panel
with open(output_sequences, "wt") as outfile:
    outfile.write("Chromosome,Position,StacksLocusID,Allele1,Allele2,NumSamples,MAF,Complexity,NumFlankingSNPs,SumFlankingMAFs,Sequence,ModifiedSequence\n")

    with open(input_filtered_info, "rt") as infile:

        for line in infile:
            chrom, pos, locus, allele1, allele2, _, _, infos = line.strip().split()
            pos = int(pos)
            locus_num, locus_pos, locus_sens = locus.split(":")
            num_samples, maf = infos.split(";")
            num_samples = num_samples.split("=")[1]
            maf = maf.split("=")[1]

            # Extract flanking sequence
            s = chromosomes[chrom][pos-flank_size-1: pos+flank_size]
            s_orig = s[:]

            # Compute complexity
            complexity = compression_ratio(s, minimum_ratio, maximum_ratio)

            # Find flanking SNPs
            flanking = []
            for i in (-1, 0, 1):
                flanking += all_snps[chrom][(pos // flank_size) + i]

            flanking = [x for x in flanking if abs(x[1] - pos) <= flank_size and abs(x[1] - pos) != 0]

            # Annotate flanking SNPs first since annotating the central SNP in the
            # [A/T] format will change the length of the string
            for f in flanking:
                a1, a2 = f[5: 7]
                pos_diff = f[1] - pos
                s = s[:flank_size+pos_diff] + "N" + s[flank_size+pos_diff+1:]

            # Then annotate central SNP
            s = s[:flank_size] + f"[{allele1}/{allele2}]" + s[flank_size+1:]

            # Report important infos
            outfile.write(",".join([
                chrom,
                str(pos),
                locus,
                allele1,
                allele2,
                num_samples,
                maf,
                str(round(complexity, 4)),
                str(len(flanking)),
                str(sum([float(x[8]) for x in flanking])),
                s_orig,
                s]) + "\n")
