#!/bin/bash
# Launch process_radtags on all the lanes

### Global variable
INFO_FILES="01-info_files"
TRIM_LENGTH=$1 # Length to cut reads after process_radtags
ENZYME=$2 # Name of the enzyme (run 'process_radtags' without options for list)

cat $INFO_FILES/lane_info.txt |
while read f
do
    grep -vE "^$" $INFO_FILES/sample_information.csv  | \
        grep "$f" | \
        grep -v "Barcode" | \
        cut -f 2 > $INFO_FILES/barcodes.txt

    # Extract the barcodes of all different lengths at once, STACKS 1.23+
    mkdir 03-samples/$f 2> /dev/null
    process_radtags \
        -i gzfastq \
        -f 02-raw/$f".fastq.gz" \
        -o 03-samples/$f \
        -b $INFO_FILES/barcodes.txt \
        -c -q -r -t $TRIM_LENGTH \
        --barcode_dist 2 \
        -E phred33 \
        -e $ENZYME

    ## Prepare bacode_lengths.txt from barcodes.txt 
    #perl -ne 'chomp; print length()."\n"' $INFO_FILES/barcodes.txt | \
    #    sort -un > $INFO_FILES/barcode_lengths.txt

    ## Create barcode files, eg: barcodes_4.txt, barcodes_5.txt ...
    #cat $INFO_FILES/barcode_lengths.txt | \
    #while read b
    #do
    #    perl -sane 'chomp; if (length eq $b) {print $_."\n"}'\
    #        -- -b=$b $INFO_FILES/barcodes.txt | \
    #    sort -u > $INFO_FILES/barcodes_$b".txt"
    #done

    ## Extract reads for each individuals using barcodes
    #cat $INFO_FILES/barcode_lengths.txt | \
    #while read b
    #do
    #    mkdir 03-samples/$f
    #    process_radtags \
    #        -i gzfastq \
    #        -f 02-raw/$f".fastq.gz" \
    #        -o 03-samples/$f \
    #        -b $INFO_FILES/barcodes_$b".txt" \
    #        -c -q -r -t $TRIM_LENGTH \
    #        --barcode_dist 2 \
    #        -E phred33 \
    #        -e $ENZYME

    #    mv 03-samples/$f/process_radtags.log 03-samples/$f/process_radtags_b$b".log" 2> /dev/null
    #    #--filter_illumina 3 \
    #done

#process_radtags [-f in_file | -p in_dir [-P] | -1 pair_1 -2 pair_2] -b barcode_file -o out_dir -e enz [-c] [-q] [-r] [-t len] [-D] [-w size] [-s lim] [-h]
#    f: path to the input file if processing single-end sequences.
#    i: input file type, either 'bustard' for the Illumina BUSTARD output files, 'fastq', or 'gzfastq' for gzipped Fastq (default 'fa stq').
#    p: path to a directory of files.
#    P: files contained within directory specified by '-p' are paired.
#    1: first input file in a set of paired-end sequences.
#    2: second input file in a set of paired-end sequences.
#    o: path to output the processed files.
#    y: output type, either 'fastq' or 'fasta' (default fastq).
#    b: path to a file containing barcodes for this run.
#    c: clean data, remove any read with an uncalled base.
#    q: discard reads with low quality scores.
#    r: rescue barcodes and RAD-Tags.
#    t: truncate final read length to this value.
#    E: specify how quality scores are encoded, 'phred33' (Illumina 1.8+, Sanger, default) or 'phred64' (Illumina 1.3 - 1.5).
#    D: capture discarded reads to a file.
#    w: set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).
#    s: set the score limit. If the average score within the sliding window drops below this value, the read is discarded

    #rm $INFO_FILES/barcode_lengths.txt 2> /dev/null
    #rm $INFO_FILES/barcodes*.txt 2> /dev/null
done

