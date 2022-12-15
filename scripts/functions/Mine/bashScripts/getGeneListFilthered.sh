                os.system("grep -v ':' rawFile.txt > loading.txt")
                os.system("grep '\.' loading.txt > rawFile.txt")
                os.system("grep -v '\[' rawFile.txt > loading.txt")
                os.system("sed 's/^.*\ //' loading.txt >> GeneList.txt")