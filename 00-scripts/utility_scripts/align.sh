for i in $(ls -1 04-all_samples/*.fq); do name=$(basename $i);bwa aln -n 3 -k 1 -t 6 ./01-info_files/Gasterosteus_aculeatus.fa $i | bwa samse -r "@RG\tID:'$name'\tSM:'$name'\tPL:Illumina" ./01-info_files/Gasterosteus_aculeatus.fa - $i > ./04-all_samples/$name.sam;done
