#!/bin/bash

# First split sample list to align into different files with:
# cd 04-all_samples
# ls -1 *.fq.gz > ../wanted_all_samples
# cd ..
# mkdir samples_split
# split -a 4 -l 41 -d wanted_all_samples samples_split/samples_split.

# Global variables
GENOMEFOLDER="08-genome"
GENOME="GCF_002021735.1_Okis_V1_genomic.fna"
DATAFOLDER="04-all_samples"
NCPU=$1
SAMPLE_FILE="$2"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi


# Index genome if not alread done
# bwa index -p $GENOMEFOLDER/$GENOME $GENOMEFOLDER/$GENOME.fasta

cat "$SAMPLE_FILE" |
while read file
do
    # Name of uncompressed file
    echo "Aligning file $file"

    name=$(basename $file)
    ID="@RG\tID:ind\tSM:ind\tPL:IonProton"

    # Align reads 1 step
    bwa mem -t "$NCPU" -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R $ID \
        $GENOMEFOLDER/$GENOME $DATAFOLDER/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > $DATAFOLDER/"${name%.fq.gz}".bam
done
