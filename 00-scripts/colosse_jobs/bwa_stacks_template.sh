#!/bin/bash
#PBS -A ihv-653-ab
#PBS -N bwa__LIST__
#PBS -o bwa__LIST__.out
#PBS -e bwa__LIST__.err
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
source /clumeq/bin/enable_cc_cvmfs
module load samtools/1.5
module load bwa/0.7.15

# Global variables
DATAFOLDER="04-all_samples"
GENOMEFOLDER="/rap/ihv-653-ab/eric_normandeau/01_projects/epic4/04_stacks_workflow_2017-06-07_2850_samples_okis_v1/08-genome/"
GENOME="GCF_002021735.1_Okis_V1_genomic.fna"

for file in $(cat $list)
do
    # Name of uncompressed file
    echo "Aligning file $file"

    name=$(basename $file)
    ID="@RG\tID:ind\tSM:ind\tPL:IonProton"

    # Align reads 1 step
    bwa mem -t 8 -k 19 -c 500 -O 0,0 -E 2,2 -T 0 \
        -R $ID \
        $GENOMEFOLDER/$GENOME $DATAFOLDER/"$name" 2> /dev/null |
        samtools view -Sb -q 1 -F 4 -F 256 -F 2048 \
        - > $DATAFOLDER/"${name%.fq.gz}".bam
done
