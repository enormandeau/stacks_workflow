#!/bin/bash
for i in `ls -1 *.fa`
do
    echo "Treating file $i"
    ../remove_mid.py $i ../barcodes.txt $i.no_barcode
done

