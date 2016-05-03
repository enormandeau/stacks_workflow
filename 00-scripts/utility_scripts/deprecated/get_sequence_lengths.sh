#!/bin/bash

for i in `ls -1 *.fa.trimmed`; do echo -ne "$i:\t"; head -20000 $i | grep -v ">" | perl -ne 'chomp; print(length."\n")' | sort | uniq -c | sort -n | tail -4; done

