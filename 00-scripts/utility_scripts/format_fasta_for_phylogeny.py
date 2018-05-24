#!/usr/bin/env python2
"""Extract wanted loci from batch_1.fa output of population
and prepare for phylogenetics analyses

- Join alleles using IUPAC degenerated nucleotides (R, Y, S, W...)
- Eliminate loci with individuals that have more than 2 alleles?
- Join all sequences to have only one long sequence per individual

Usage:
    ./00-scripts/utility_scripts/format_fasta_for_phylogeny.py input_fasta wanted_loci output_fasta
"""

# Modules
from collections import defaultdict
import sys

# Classes
class Fasta(object):
    """Fasta object with name and sequence
    """
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

        info = name.split("_")
        self.locus = int(info[1])
        self.allele = int(info[7].split(" ")[0])

        info = name.split(" ")
        self.sample = info[1].replace("[", "").replace("]", "")

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

    def __repr__(self):
        return "\n".join([self.name,
                        self.sequence,
                        str(self.locus),
                        str(self.allele),
                        self.sample])

# Functions
def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator
    """
    with open(input_file) as f:
        sequence = ""
        name = ""
        begun = False
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if begun:
                    yield Fasta(name, sequence)
                name = line.replace(">", "")
                sequence = ""
                begun = True
            else:
                sequence += line
        yield Fasta(name, sequence)

def concensus(s1, s2):
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

def write_output_fasta(output_fasta, seq_dict, sample_set):
    sample_sequences = defaultdict(list)
    with open(output_fasta, "w") as outfile:
        for locus in seq_dict:
            allele_count = defaultdict(int)
            for sample in sorted(sample_set):
                try:
                    sample_sequences[sample].append(seq_dict[locus][sample])
                    seq = seq_dict[locus][sample]
                    allele_count[seq] += 1
                except:
                    sample_sequences[sample].append("NA")

            # Find most frequen allele
            most_frequent_count = max(allele_count.values())
            most_frequent = [x for x in allele_count.items() if x[1] == most_frequent_count][0][0]

            # Replace missing data by most frequent allele
            for sample in sample_set:
                if sample_sequences[sample][-1] == "NA":
                    sample_sequences[sample][-1] = most_frequent

        # Concatenate all sequences for each sample and write to output file
        for sample in sample_sequences:
            seq = "".join(sample_sequences[sample])
            outfile.write(">" + sample + "\n" + seq + "\n")

# Parse user input
try:
    input_fasta = sys.argv[1]
    wanted_loci = sys.argv[2]
    output_fasta = sys.argv[3]
except:
    print __doc__
    sys.exit(1)

# Program
seq_dict = dict()
sequences = fasta_iterator(input_fasta)

# Get wanted loci numbers
wanted = set()
with open(wanted_loci) as infile:
    for line in infile:
        wanted.add(int(line.strip()))

# Get wanted sequences
sample_set = set()
for seq in sequences:
    if seq.locus in wanted:
        sample_set.add(seq.sample)
        if seq.locus in seq_dict and seq.allele < 2:
            if seq.sample in seq_dict[seq.locus]:
                seq_dict[seq.locus][seq.sample] = concensus(seq_dict[seq.locus][seq.sample], seq.sequence)
            else:
                seq_dict[seq.locus][seq.sample] = seq.sequence
        else:
            seq_dict[seq.locus] = dict()
            seq_dict[seq.locus][seq.sample] = seq.sequence

# Write output
write_output_fasta(output_fasta, seq_dict, sample_set)
