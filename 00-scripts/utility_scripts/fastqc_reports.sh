#!/bin/bash

# Global variables
INPUTFOLDER=$1
DESCRIPTION=$(basename $INPUTFOLDER)
OUTPUTFOLDER="90-fastqc_$DESCRIPTION"

# Test if folder exists
if ! [[ -d "$INPUTFOLDER" ]]
then
    echo "Usage:"
    echo "    ./00-scripts/utility_scripts/fastqc_reports.sh INPUTFOLDER"
    exit 1
fi

# Create fastqc_output_folder
mkdir "$OUTPUTFOLDER" 2> /dev/null

# iterate through fastq.gz or fq.gz files
for file in "$INPUTFOLDER"/*.f*q.gz
do
    tempfile="$OUTPUTFOLDER"/"$(basename $file)"
    gunzip -c "$file" | head -100000 | gzip -c > "$tempfile"
    fastqc "$tempfile"
done

# Create multiqc report
if hash multiqc 2>/dev/null; then
    echo "running multiqc to create quality report"
    multiqc -f -i "$DESCRIPTION" -n fastqc_reports_"$DESCRIPTION".html "$OUTPUTFOLDER"
else
    echo "install multiqc to get a quality report"
fi
