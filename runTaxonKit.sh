#!/bin/bash

if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

taxonkit list --ids 9443 -nr --indent "" > Output.txt
grep -r "species" ./Output.txt > 1.txt
grep -v "subspecies" ./1.txt > Output.txt
sed 's/\s.*$//' Output.txt  > SpeciesIdsFromRank.txt
rm 1.txt
