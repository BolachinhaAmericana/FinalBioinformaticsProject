grep -v ':' rawFile.txt > loading1.txt # Remove Lines with :
grep '\.' loading1.txt > rawFile1.txt # filthers only lines with '. '
grep -v '\[' rawFile1.txt > loading2.txt # sometimes there where lines with more then just the gene name. So we just remove everything with '['
sed 's/^.*\ //' loading2.txt >> GeneList.txt # removes everything behind (including) '. '