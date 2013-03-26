#!/bin/bash

cat 01-info_files/lane_info.txt | while read f; do for s in $(ls -1 03-samples/$f/sample*); do sample=$(basename $s); mv 03-samples/$f/$sample 03-samples/$f/$f"_"$sample; done; done

cat 01-info_files/lane_info.txt | while read f; do for s in $(ls -1 03-samples/$f/*.fq); do sample=$(basename $s); cp -l 03-samples/$f/$sample 04-all_samples; done; done

