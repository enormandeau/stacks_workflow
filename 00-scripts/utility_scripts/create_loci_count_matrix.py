#!/usr/bin/env python
"""Create matrix of number of reads per catalog locus
from *.matches.tsv.gz files in 05_stacks or in 05_stacks_rx

Usage:
    ./00-scripts/utility_scripts/create_matches_matrix.py input_folder presence_threshold output_file
"""

# Modules
from collections import defaultdict
import copy
import gzip
import sys
import os

# Classes

# Functions

# Main
if __name__ == '__main__':
    try:
        input_folder = sys.argv[1]
        presence_threshold = float(sys.argv[2])
        output_file = sys.argv[3]
    except:
        print __doc__
        sys.exit(1)

    matches_files = sorted([x for x in os.listdir(input_folder) if x.endswith(".matches.tsv.gz")])

    # TESTING
    #matches_files = matches_files[0:4]

    print "There are:", len(matches_files), "matches files to treat"

    matches_dict = defaultdict(lambda: defaultdict(int))
    ind = 1
    num_ind = len(matches_files)
    threshold = presence_threshold * float(num_ind)
    samples = []

    for mfile in matches_files:
        sample_name = mfile.replace(".matches.tsv.gz", "")
        samples.append(sample_name)
        print "  Treating individual {}/{}: {}".format(ind, num_ind, sample_name)
        ind += 1

        m = gzip.open(os.path.join(input_folder, mfile))

        # Getting rid of first line
        m.readline()

        # Populating matches_dict
        for line in m:
            l = line.strip().split()
            locus, count = l[2], int(l[6])

            matches_dict[locus][sample_name] += count

        m.close()

    # Trim down the matches_dict with threshold
    print "Keeping only loci with at least {} individuals".format(threshold)
    to_remove = []

    for i in matches_dict:
        if len(matches_dict[i]) < threshold:
            to_remove.append(i)

    trimmed_dict = copy.deepcopy(matches_dict)
    for i in to_remove:
        trimmed_dict.pop(i)

    print "  Reduced number of loci by {:2.2f}% ({} to {} loci)".format(100.0 - 100.0 * len(trimmed_dict) / len(matches_dict), len(matches_dict), len(trimmed_dict))

    # Normalize matrix per individual depth
    sample_depth = dict()
    for sample in samples:
        total = 0
        for locus in trimmed_dict:
            total += trimmed_dict[locus][sample]

        sample_depth[sample] = float(total)

    for sample in samples:
        for locus in trimmed_dict:
            trimmed_dict[locus][sample] = 1000000. * trimmed_dict[locus][sample] / sample_depth[sample]

    # Writing matrix of counts to file
    outfile = open(output_file, "w")

    outfile.write("\t".join(samples) + "\n")

    for locus in trimmed_dict:
        counts = []
        for sample in samples:
            counts.append(str(trimmed_dict[locus][sample]))

        outfile.write("\t".join(counts) + "\n")

    outfile.close()

