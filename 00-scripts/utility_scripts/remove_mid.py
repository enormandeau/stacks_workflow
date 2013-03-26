#!/usr/bin/python
"""
Remove MID sequences found at the beginning of each sequence

Usage:
    remove_mid.py  sequence_file  mid_file  output_file

The sequence file contains sequences in Fasta format

The MID file contains one MID sequence per line and nothing else. Be sure to put
the longuest MIDs first in the file so as to avoid the very possible case where
shorter MID sequences are found within longer ones.
"""

# Importing libraries
import sys
try:
    from Bio import SeqIO
except:
    print """This program requires the Biopython library
    In Ubutu or Debian, install with:
    sudo apt-get install python-biopython"""
    sys.exit(0)

# Defining function
def get_mids(mid_file):
    """Read file containing MIDs and put them in a list
    """
    mids = []
    with open(mid_file) as f:
        for line in f:
            l = line.strip()
            if l != "":
                mids.append(l)
    return mids

def remove_mids(sequence_file, mids):
    """Treat sequence file to remove MIDs
    """
    fasta_sequences = SeqIO.parse(open(sequence_file),'fasta')
    with open(output_file, "w") as out_f:
        for s in fasta_sequences:
            name = s.id
            sequence = s.seq.tostring()
            for m in mids:
                if sequence.startswith(m):
                    length = len(m)
                    sequence = sequence[length:]
                    out_f.write(">" + name + "\n" + sequence + "\n")
                    break

# Parsing command line options
try:
    sequence_file = sys.argv[1]
    mid_file = sys.argv[2]
    output_file = sys.argv[3]
except:
    print __doc__
    sys.exit(1)

try:
    with open(sequence_file) as f:
        pass
except:
    print "The sequence file was not found"
    sys.exit(1)

try:
    with open(mid_file) as f:
        pass
except:
    print "The MID file was not found"
    sys.exit(1)

# Removing mids
mids = get_mids(mid_file)
remove_mids(sequence_file, mids)

