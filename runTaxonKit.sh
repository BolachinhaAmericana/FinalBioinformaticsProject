taxonkit list --ids 9443 -nr --indent "" > Output.txt
grep -r "species" ./Output.txt > 1.txt
grep -v "subspecies" ./1.txt > Output.txt
sed 's/\s.*$//' Output.txt  > SpeciesIdsFromRank.txt
rm 1.txt
rm Output.txt
