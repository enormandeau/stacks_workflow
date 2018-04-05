#!/bin/bash

# Global variables
GENOMEFOLDER="08-genome"
GENOME="kelpfly_genome_AA_pacbio_100kb_plus.fasta"
DATAFOLDER="04-all_samples"
NCPU=$1

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=4
fi

# Index genome if not alread done
# bwa index -p $GENOMEFOLDER/$GENOME $GENOMEFOLDER/$GENOME.fasta

module load bwa
module load samtools

for file in $(ls -1 $DATAFOLDER/*_R1_*.fastq.gz)
do
    # Name of uncompressed file
    file2=$(echo "$file" | perl -pe 's/_R1_/_R2_/')
    echo "Aligning file $file $file2" 

    name=$(basename $file)
    name2=$(basename $file2)
    ID="@RG\tID:ind\tSM:ind\tPL:IonProton"

    # Align reads 1 step
    bwa mem -t "$NCPU" -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R $ID \
        $GENOMEFOLDER/$GENOME $DATAFOLDER/"$name" $DATAFOLDER/"$name2" 2> /dev/null |
        samtools view -Sb -q 20 \
            -f 83 -f 163 -f 99 -f 147 \
        - > $DATAFOLDER/"${name%.fq.gz}".bam

    # Sort and index
    samtools sort $DATAFOLDER/"${name%.fq.gz}".bam > $DATAFOLDER/"${name%.fq.gz}".sorted.bam
    samtools index $DATAFOLDER/"${name%.fq.gz}".sorted.bam
done
