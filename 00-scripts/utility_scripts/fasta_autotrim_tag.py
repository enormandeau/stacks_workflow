#!/usr/bin/python
"""Automatically trim the remaining enzyme site and tag nucleotides from fasta file

Usage:
    %program  inputfile  outputfile  [minlength]  [maxlength]  [minpercent]
"""

import sys
import itertools
from collections import defaultdict
from Bio import SeqIO

try:
    from Bio import SeqIO
except:
    print "This program requires the Biopython library"
    sys.exit(0)

try:
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
except:
    print __doc__
    print "You must give both and an input and an output file"
    sys.exit(1)

try:
    minlength=int(sys.argv[3])
except:
    minlength = 1

try:
    maxlength=int(sys.argv[4])
except:
    maxlength = 20

try:
    minpercent=int(sys.argv[5])
except:
    minpercent = 75

#print "Parameters:"
#print "minlength:", minlength
#print "maxlength:", maxlength
#print "minpercent:", minpercent, "\n"

class Indexable(object):
    def __init__(self,it):
        self.it=it
        self.already_computed=[]
    def __iter__(self):
        for elt in self.it:
            self.already_computed.append(elt)
            yield elt
    def __getitem__(self,index):
        try:
            max_idx=index.stop
        except AttributeError:
            max_idx=index
        n=max_idx-len(self.already_computed)+1
        if n>0:
            self.already_computed.extend(itertools.islice(self.it,n))
        return self.already_computed[index]

print inputfile, outputfile, minlength, maxlength, minpercent

handle = open(inputfile, 'rU')
sequences = Indexable(SeqIO.parse(handle, 'fasta'))[0:10000]
seqs = [s.seq.tostring() for s in sequences]
handle.close()

tag_length = -99

for i in xrange(minlength, maxlength + 1):
    d = defaultdict(int)
    #print i
    short = [s[0:i] for s in seqs]
    for s in short:
        d[s] += 1
    frequencies = d.values()
    #print frequencies
    percent = 100.0 * max(frequencies) / sum(frequencies)
    if percent <= minpercent:
        tag_length = i - 1
        break

if tag_length > 0:
    print "   >>> Tag length is:", tag_length
    with open(outputfile, "w") as f:
        handle = open(inputfile, 'rU')
        for s in SeqIO.parse(handle, 'fasta'):
            f.write(">" + s.id + "\n")
            f.write(s.seq.tostring()[tag_length:] + "\n")
    #print "   >>> Now trimming file:", inputfile
else:
    print "Did not find tag length. Try higher values for maxlength."
    sys.exit(1)






































