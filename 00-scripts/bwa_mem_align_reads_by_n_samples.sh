#!/bin/bash

# First split sample list to align into different files with:
# cd 04-all_samples
# ls -1 *.fq.gz > ../all_samples_for_alignment.txt
# cd ..
# mkdir samples_split
# split -a 4 -l 1 -d all_samples_for_alignment.txt samples_split/samples_split.

## With GNU Parallel
# ls -1 samples_split/* | parallel -k -j 10 ./00-scripts/bwa_mem_align_reads_by_n_samples.sh 4 {}

## With GNU Parallel on slurm
# ls -1 samples_split/* | parallel -k -j 10 srun -c 4 --mem 20G -p large --time 21-00:00 -J bwaMem -o 10-log_files/bwaMEMsplit_%j.log ./00-scripts/bwa_mem_align_reads_by_n_samples.sh 4 {} \; sleep 0.1 &

## With srun on a single file
# srun -c 4 --mem 20G -p large --time 21-00:00 -J bwaMem -o 10-log_files/bwaMEMsplit_%j.log ./00-scripts/bwa_mem_align_reads_by_n_samples.sh 4 <SAMPLE_FILE>

# Global variables
GENOMEFOLDER="08-genome"
GENOME="genome.fasta"
DATAFOLDER="04-all_samples"
NCPU="$1"
SAMPLE_FILE="$2"

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

cat "$SAMPLE_FILE" |
while read file
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
