#!/bin/bash
# Combine read files from two folders where the two folders contain files with
# the same names

# Get user input
dir1=$1
dir2=$2
dir3=$3

# Create output directory if it does not exist
mkdir $dir3

# Combine the files
for i in d1/*.txt; do file=$(basename $i); cat $dir1/$file $dir2/$file > $dir3/$file; done

