#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N gsnap__LIST__
#PBS -o gsnap__LIST__.out
#PBS -e gsnap__LIST__.err
#PBS -l walltime=24:00:00
#PBS -M your.email.address@here
#PBS -m bea
#PBS -l nodes=1:ppn=8
#PBS -r n

# Global variables
list=__LIST__

# Move to job submission directory
cd $PBS_O_WORKDIR

# Load gmap module
module load apps/gmap/2015-12-31.v9

#prepare the genome if not already done
#gmap_build --dir=/rap/ihv-653-ab/00_ressources/01_genomes/Oniloticus/ /rap/ihv-653-ab/00_ressources/01_genomes/Oniloticus/oni_ref_Orenil1.1.fa -d gmap_oniloticuv1.1

# Global variables
DATAFOLDER="04-all_samples"
GENOMEFOLDER="/rap/ihv-653-ab/00_ressources/01_genomes/Omykiss/"
GENOME="gmap_omykiss"
PWD="__PWD__"

for file in $(cat $list)
do
    # Align reads
    echo "Aligning $file"
    gsnap --gunzip -t 8 -A sam -m 1 -i 2 --min-coverage=0.90 \
        --dir="$GENOMEFOLDER" -d "$GENOME" \
        -o "${file%.fq.gz}".sam \
        "$file"

    # Create bam file
    echo "Creating bam for $file"
    samtools view -Sb -q 1 -F 4 -F 1797 \
        "${file%.fq.gz}".sam > "${file%.fq.gz}".bam

    # Clean up
    echo "Removing $file"
    rm $DATAFOLDER/"${file%.fq.gz}".sam
done
