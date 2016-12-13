#!/bin/bash

# Global variables
GENOMEFOLDER="08-genome"
GENOME="genome.fasta"
DATAFOLDER="04-all_samples"
NCPU=$1

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi

# Index genome if not alread done
# bwa index -p $GENOMEFOLDER/$GENOME $GENOMEFOLDER/$GENOME.fasta

for file in $(ls -1 $DATAFOLDER/*.fq.gz)
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
