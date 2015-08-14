#!/usr/bin/env python
"""Combine distribution graphs of two vcf files to compare them

Usage:
    ./combine_distribution_graphs.py folder_1 folder_2 output_folder

Where the folders_N are the output folders of 00-scripts/05_filter_vcf.py used
with the -g option and the output_folder is the desired output folder
"""

# Modules
import Image
import sys
import os

# Function
def combine_images(image_1, image_2, image_out, xdim=800, ydim=800):
    """Combine two images vertically
    """
    i1 = Image.open(image_1)
    i2 = Image.open(image_2)
    new_image = Image.new("RGB", (800, 800))

    new_image.paste(i1, (0, 0))
    new_image.paste(i2, (0, 400))
    new_image.save(image_out)

def error(message):
    print("Error: {}".format(message))

if __name__ == '__main__':
    # Parse user input
    try:
        folder_1 = sys.argv[1]
        folder_2 = sys.argv[2]
        output_folder = sys.argv[3]
    except:
        print __doc__
        sys.exit(1)

    # Confirm the folders exist
    if not os.path.exists(folder_1):
        error("Could not find folder {}".format(folder_1))
        sys.exit(1)

    if not os.path.exists(folder_2):
        error("Could not find folder {}".format(folder_2))
        sys.exit(1)

    # Create the output folder and subdirectories
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    subdirectory_1 = os.path.join(output_folder, "global")
    subdirectory_2 = os.path.join(output_folder, "populations")

    if not os.path.exists(subdirectory_1):
        os.makedirs(subdirectory_1)

    if not os.path.exists(subdirectory_2):
        os.makedirs(subdirectory_2)

    # For all images in folder_1, merge with equivalent in folder_2
    wanted = dict()
    wanted["global"] = os.listdir(os.path.join(folder_1, "global"))
    wanted["populations"] = os.listdir(os.path.join(folder_1, "populations"))

    for w in wanted:
        for f in wanted[w]:
            i1 = os.path.join(folder_1, w, f)
            i2 = os.path.join(folder_2, w, f)
            iout = os.path.join(output_folder, w, f)
            combine_images(i1, i2, iout)

