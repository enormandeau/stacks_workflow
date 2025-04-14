#!/usr/bin/env python3
"""Extract wanted loci from a filtered VCF and sequences from populations.loci.fa and
prepare for phylogenetics analyses

- Join alleles using IUPAC degenerated nucleotides (R, Y, S, W...)
- Join all sequences to have only one long sequence per individual
- Input and output files can be gzip compressed

Usage:
    ./00-scripts/utility_scripts/format_fasta_for_phylogeny_stacks2.py input_vcf input_fasta input_genome output_fasta
"""

# Modules
from collections import defaultdict
import gzip
import sys

# Classes
class Fasta(object):
    """Fasta object with name and sequence
    """
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

def consensus(s1, s2):
    result = []
    for i in range(len(s1)):
        n1 = s1[i]
        n2 = s2[i]
        if n1 == n2:
            result.append(n1)
        else:
            result.append(degenerated(n1, n2))

    result = "".join(result)
    return result

def degenerated(n1, n2):
    iupac_table = {
            "AC": "M",
            "AG": "R",
            "AT": "W",
            "CG": "S",
            "CT": "Y",
            "GT": "K"}

    return iupac_table["".join(sorted([n1, n2]))]

# Parse user input
try:
    input_vcf = sys.argv[1]
    input_fasta = sys.argv[2]
    input_genome = sys.argv[3]
    output_fasta = sys.argv[4]
except:
    print(__doc__)
    sys.exit(1)

# Store loci infos and sequences from populations.loci.fa
locus_dict = dict()
sequences = fasta_iterator(input_fasta)

for s in sequences:
    locus_id, scaf, pos, sense = s.name.replace("CLocus_", "").replace("[", "").replace("]", "").replace(",", "").split(" ")
    s.name = " ".join([locus_id, scaf, pos])
    pos = int(pos)
    locus_dict[locus_id] = s

# Store SNP infos from VCF
snp_dict = defaultdict(list)

with myopen(input_vcf, "rt") as infile:
    for line in infile:
        if line.startswith("##"):
            continue
        elif line.startswith("#CHROM"):
            sample_names = line.strip().split("\t")[9: ]
        else:
            l = line.strip().split("\t")
            locus_id, pos, sense = l[2].split(":")
            pos = int(pos)
            alleles = "".join(l[3: 5])
            genotypes = [x.split(":")[0].replace("/", "") for x in l[9: ]]

            snp_dict[locus_id].append((pos, sense, alleles, genotypes))

# Extract regions from genome
scaffolds_dict = dict()
genome = fasta_iterator(input_genome)

for s in genome:
    shortname = s.name.split(" ")[0]
    scaffolds_dict[shortname] = s.sequence.upper()

# For each locus_id in snp_dict, get sequence from locus_dict
sequences_per_sample = defaultdict(list)

for l in snp_dict:
    locus = locus_dict[l]
    locus_id, scaf, pos = locus.name.split(" ")
    pos = int(pos)
    snps = snp_dict[l]

    # Get sequence from genome
    if snps[0][1] == "+":
        seq = scaffolds_dict[scaf][pos: pos+130]

    elif snps[0][1] == "-":
        seq = scaffolds_dict[scaf][pos-130: pos]

    if len(seq) < 130:
        continue

    # Change position in sequence
    if snps[0][1] == '+':
        #print(seq[snps[0][0]-2] in snps[0][2])
        seq_pos = snps[0][0] - 2

    elif snps[0][1] == '-':
        #print(seq[-snps[0][0]] in snps[0][2])
        seq_pos = snps[0][0]

    for i, sample in enumerate(sample_names):
        seq = []
        for snp in snps:
            alleles = snp[2]
            genotypes = snp[3]
            n1 = alleles[int(genotypes[i][0])]
            n2 = alleles[int(genotypes[i][1])]
            cons = consensus(n1, n2)

            # Modify sequence
            #seq = seq[:seq_pos] + cons + seq[seq_pos+1:]
            seq.append(cons)

        sequences_per_sample[sample].append("".join(seq))

# Write all the sequences
with myopen(output_fasta, "wt") as outfile:
    for sample in sequences_per_sample:
        outfile.write(">" + sample + "\n")
        outfile.write("".join(sequences_per_sample[sample]) + "\n")
