#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N gsnap__LIST__
#PBS -o gsnap__LIST__.out
#PBS -e gsnap__LIST__.err
#PBS -l walltime=20:00:00
#PBS -M YOUREMAIL
####PBS -m ea
#PBS -l nodes=1:ppn=8
#PBS -r n

# Global variables
list=__LIST__

# Move to job submission directory
cd $PBS_O_WORKDIR

# Load gmap module
module load apps/gmap/2015-12-31.v9
module load compilers/gcc/4.8.5  apps/mugqic_pipeline/2.1.1
module load mugqic/samtools/1.2

#prepare the genome if not already done
#gmap_build --dir=/rap/ihv-653-ab/00_ressources/01_genomes/Oniloticus/ /rap/ihv-653-ab/00_ressources/01_genomes/Oniloticus/oni_ref_Orenil1.1.fa -d gmap_oniloticuv1.1

# Global variables
DATAFOLDER="04-all_samples"
GENOMEFOLDER="/rap/ihv-653-ab/00_ressources/01_genomes/Omykiss/"
GENOME="gmap_omykiss"

for file in $(cat $list)
do
    # Align reads
    echo "Aligning $file"
    gsnap --gunzip -t 8 -A sam -m 5 -i 2 --min-coverage=0.90 \
        --dir="$GENOMEFOLDER" -d "$GENOME" \
        --read-group-id="${file%.fq.gz}" \
        -o "${file%.fq.gz}".sam \
        "$file"

    # Create bam file
    echo "Creating bam for $file"
    samtools view -Sb -q 1 -F 4 -F 256 \
        "${file%.fq.gz}".sam > "${file%.fq.gz}".bam

    # Clean up
    echo "Removing $file"
    rm "${file%.fq.gz}".sam
done
