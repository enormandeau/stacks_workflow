#!/bin/bash
# Replace string 'your.email@service.com' by content of $YOUREMAILADDRESS
# Set this last variable in your ~/.bashrc file

for script in 00-scripts/katak_jobs/*.sh
do
    grep "your.email@service.com" "$script"
    perl -i -sape 's/YOUREMAIL/$e/' -- -e="$YOUREMAIL" "$script"
done
