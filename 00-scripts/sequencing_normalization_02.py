#!/usr/bin/env python
"""Calculate new volumes needed to normalize the sequencing depth per sample
and then create a DNA plate template with these volumes for a new library.

Usage:
    ./00-scripts/sequencing_normalization_02.py targetNumReads minimumReads totalVolume

Where:
    targetNumReads = desired total number of sequenced reads [int]
"""

# Modules
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os

# Classes

# Functions

# Parsing user input
try:
    targetNumReads = int(sys.argv[1])
    minimumReads = int(sys.argv[2])
    totalVolume = float(sys.argv[3])
except:
    print __doc__
    sys.exit(1)

# Main

# Global variables
folder="03-samples"

chips = [x.replace(".infos", "") for x in os.listdir(folder) if x.endswith(".infos")]

# Iterate over the chips
for chip in chips:

    # File names
    info_file = os.path.join(folder, chip + ".infos")
    numseq_file = os.path.join(folder, chip + ".numseq")

    # Load data into pandas data frame
    info_df = pd.read_csv(info_file, sep="\t", header=None)
    info_df = info_df.ix[:,0:5]
    info_df.columns = ["Chip", "Barcode", "Population", "Individual", "PopID", "Well"]

    numseq_df = pd.read_csv(numseq_file, sep=" ", header=None)
    numseq_df.columns = ["Barcode", "NumReads"]

    # Merge infos and number of reads into one dataframe
    data = pd.merge(info_df, numseq_df, on="Barcode")

    # Computations
    data["Missing"] = targetNumReads - data["NumReads"]
    data.loc[data["NumReads"] > targetNumReads, "Missing"] = 0
    data.loc[data["NumReads"] < minimumReads, "Missing"] = 0

    data["Correction"] = data["Missing"].astype(float) / data["NumReads"]

    data["Volume"] = data["Correction"] / sum(data["Correction"]) * totalVolume
    data.loc[data["Volume"] < 1, "Volume"] = 1.0

    # Create output csv file
    rows = list("ABCDEFGH")
    columns = range(1, 13)

    plate = pd.DataFrame(np.zeros([8, 12]))
    plate.columns = columns
    plate.index = rows

    for well in data["Well"]:
        row = well[0]
        column = int(well[1:])
        volume = float(data.loc[data["Well"] == well, "Volume"])
        volume = "{0:.1f}".format(volume)
        plate.loc[row, column] = volume

    print plate
    plate.to_csv(os.path.join(folder, chip + "_normalization.csv"), sep="\t")
