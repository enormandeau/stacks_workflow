#!/bin/bash
# Renaming the extracted sample files
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename "$0")
LOG_FOLDER="10-log_files"

cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Global variables
INFO_FILES="01-info_files"
SAMPLES_FOLDER="03-samples"
ALL_SAMPLES_FOLDER="04-all_samples"

# Renaming files
grep -vE '^#' "$INFO_FILES"/sample_information.csv |
    cut -f 1-4 |
    perl -pe 's/\.f(ast)*q\.gz//' |
    perl -pe 's/\t/\/sample_/' |
    perl -pe 's/([ACTG]+)\t/\1.fq.gz\t/' |
    awk '{print $1"\t04-all_samples/"$2"_"$3".fq.gz"}' > renaming_01.txt
    #perl -pe 's/(\w+)\t(\S+)$/04-all_samples\/\1_\2.fq.gz/' #> renaming_01.txt

cut -f 2 renaming_01.txt | sort -u > renaming_02.txt

cat renaming_02.txt |
    while read -r i
    do
        echo Treating: "$i"
        rm $ALL_SAMPLES_FOLDER/"$i" 2> /dev/null
        num_copies=$(grep "$i" renaming_01.txt | cut -f 1 | wc -l)
        echo -n "  Sample found $num_copies times: "
        if [[ $num_copies -eq 1 ]]
        then
            echo "Creating a link"
            grep "$i" renaming_01.txt |
                cut -f 1 |
                while read -r j
                do
                    echo "    Linking: $SAMPLES_FOLDER/""$j"
                    ln $SAMPLES_FOLDER/"$j" "$i"
                done
        else
            echo "Merging samples"
            grep "$i" renaming_01.txt |
                cut -f 1 |
                while read -r j
                do
                    echo "    Copying: $SAMPLES_FOLDER/""$j"
                    cat $SAMPLES_FOLDER/"$j" >> "$i"
                done
        fi
    done | tee 10-log_files/"$TIMESTAMP"_03_rename_samples_complex.log

rm renaming_01.txt renaming_02.txt 2> /dev/null

