#!/bin/bash

# Global variables
DATAFOLDER="04-all_samples"
GENOMEFOLDER="07-genome"
GENOME="sfon_genome_v1.0.fasta"

# Index genome if not alread done
# bwa index -p $GENOMEFOLDER/$GENOME $GENOMEFOLDER/$GENOME.fasta

for file in $(ls -1 $DATAFOLDER/*.fq.gz)
do
    # Name of uncompressed file
    i=${file%.gz}
    echo "Treating file $i"

    # Decompress file
    echo "  Decompressing fastq file..."
    gunzip -c $file > $i

    name=$(basename $i)
    ind=$(echo $name | cut -d "_" -f 2)
    ID="@RG\tID:ind${ind}\tSM:ind${ind}\tPL:IonProton"

    # Align reads
    # NOTES:
    # test -c 1 10 100 1000
    # test different k values
    # -L 100,100 to avoid soft clipping -> does it work?
    echo "  Aligning $name"
    bwa mem -t 16 -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R $ID \
        $GENOMEFOLDER/$GENOME $DATAFOLDER/"$name" > $DATAFOLDER/"${name%.fq}".sam \
        2> /dev/null

    # Create bam file
    echo "  Creating bam file..."
    samtools view -Sb -q 1 -F 4 -F 256 -F 1797 \
        $DATAFOLDER/"${name%.fq}".sam > $DATAFOLDER/"${name%.fq}".unsorted.bam

    # Sort and index bam file
    echo "  Sorting bam file..."
    samtools sort $DATAFOLDER/"${name%.fq}".unsorted.bam $DATAFOLDER/"${name%.fq}"

    echo "  Indexing bam file..."
    samtools index $DATAFOLDER/"${name%.fq}".bam

    # Clean up
    rm $DATAFOLDER/"${name%.fq}".sam
    rm $DATAFOLDER/"${name%.fq}".unsorted.bam
    rm $i # uncompressed file

done
