#!/bin/bash
for i in *.bam; do echo "$i $(samtools view $i | wc -l) $(gunzip -c ${i%.sorted.bam}.fq.gz |wc -l | awk '{print $1/4}')"; done | tee ../num_sequences_per_bam_and_fqgz.txt
