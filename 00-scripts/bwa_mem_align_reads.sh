#!/bin/bash

# Global variables
DATAFOLDER="04-all_samples"
GENOMEFOLDER="07-genome"
GENOME="parastichopus_parvimensis_genome_contigs"

# Index genome if not alread done
# bwa index -p $GENOMEFOLDER/$GENOME $GENOMEFOLDER/$GENOME.fasta

for file in $(ls -1 $DATAFOLDER/*.fq.gz)
do
    # Name of uncompressed file
    i=${file%.gz}
    echo "Decompression $i"

    # Decompress file
    gunzip -c $file > $i

    name=$(basename $i)
    ind=$(echo $name | cut -d "_" -f 2)
    ID="@RG\tID:ind${ind}\tSM:ind${ind}\tPL:IonProton"
    echo $name
    echo $DATAFOLDER/"$name"

    # Align reads
    bwa mem -t 16 -k 19 -c 500 -O 0,0 -E 2,2 -v 1 -T 10 -a -h 3 -a \
        -R $ID \
        $GENOMEFOLDER/$GENOME $DATAFOLDER/"$name" > $DATAFOLDER/"${name%.fq}".sam

    # Create bam file
    samtools view -Sb $DATAFOLDER/"${name%.fq}".sam > $DATAFOLDER/"${name%.fq}".unsorted.bam

    # Sort and index bam file
    samtools sort $DATAFOLDER/"${name%.fq}".unsorted.bam $DATAFOLDER/"${name%.fq}".bam
    samtools index $DATAFOLDER/"${name%.fq}".bam

    ## Clean up
    #rm $DATAFOLDER/"${name%.fq}".sam
    #rm $DATAFOLDER/"${name%.fq}".unsorted.bam
    #rm $i

done

