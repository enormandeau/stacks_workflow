#!/bin/bash
# srun -c 4 --mem 20G -p large --time 21-00:00 -J bwaMem -o 10-log_files/bwaMEM_%j.log ./00-scripts/bwa_mem_align_reads.sh 4

# Global variables
GENOMEFOLDER="08-genome"
GENOME="genome.fasta"
DATAFOLDER="04-all_samples"
NCPU="$1"

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=1
fi

# Modules
module load bwa
module load samtools

# Index genome if not alread done
# bwa index "$GENOMEFOLDER"/"${GENOME%.fasta}"

for file in $(ls -1 "$DATAFOLDER"/*.fq.gz)
do
    # Name of uncompressed file
    echo "Aligning file $file"

    name=$(basename "$file")
    ID="@RG\tID:ind\tSM:ind\tPL:IonProton"

    # Align reads 1 step
    bwa mem -t "$NCPU" -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R "$ID" \
        "$GENOMEFOLDER"/"$GENOME" "$DATAFOLDER"/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > "$DATAFOLDER"/"${name%.fq.gz}".bam

    # Samtools sort
    samtools sort --threads "$NCPU" -o "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam \
        "$DATAFOLDER"/"${name%.fq.gz}".bam

    samtools index "$DATAFOLDER"/"${name%.fq.gz}".sorted.bam

    # Cleanup
    rm "$DATAFOLDER"/"${name%.fq.gz}".bam
done
