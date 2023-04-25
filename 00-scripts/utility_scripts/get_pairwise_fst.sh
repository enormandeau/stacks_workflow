#!/bin/bash

# Inputs
VCF="$1"    # Compressed VCF
POPS="$2"   # File with one population name per line
OUT="$3"    # Output file name

# Prepare output file
rm "$OUT"
echo -e "Pop1\tPop2\tFst" > "$OUT"

# Prepare files with sample IDs
cat "$POPS" |
    while read i
    do
        echo "$i"
        grep "$i" popmap_for_admixture_plot.tsv | cut -f 1 > population_"$i".ids
    done

# Get Fst values with vcftools
for i in $(cat "$POPS")
do for j in $(cat "$POPS")
    do
        if [ $i = $j ]
            then break
        fi

        echo Treating "$i" "$j"

        vcftools \
            --gzvcf "$VCF" \
            --weir-fst-pop population_"$i".ids \
            --weir-fst-pop population_"$j".ids \
            --out "$i"-"$j"

        # Output results
        echo -e "$i\t$j\t"$(grep -v ^CHROM "$i"-"$j".weir.fst | cut -f 3 | grep -v nan | sort -n | awk '{t+=$1} END {print t/NR}' ) >> "$OUT"

    done
done

column -t "$OUT"
