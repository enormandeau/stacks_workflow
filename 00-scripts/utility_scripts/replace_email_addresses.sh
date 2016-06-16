#!/bin/bash
# Replace string 'your.email@service.com' by content of $YOUREMAILADDRESS
# Set this last variable in your ~/.bashrc file

for script in 00-scripts/katak_jobs/*.sh
do
    echo "$script"
    grep "your.email@service.com" "$script"
    echo "$YOUREMAIL"
    perl -i.bak -sape 's/YOUREMAIL/$e/' -- -e="$YOUREMAIL" "$script"
done
