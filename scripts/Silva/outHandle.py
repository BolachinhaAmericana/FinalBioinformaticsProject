from Bio import Entrez
import os
import sys

Entrez.email = "guilherme.vcc.silva@gmail.com" 

db= sys.argv[1] 
term = sys.argv[2]


search_handle = Entrez.esearch(db=db, term=term, usehistory="y", idtype="acc")
search_results = Entrez.read(search_handle)
search_handle.close()

acc_list = search_results["IdList"]
print(len(acc_list))
count = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

batch_size = 9999

for start in range(0, count, batch_size):
    end = count
    print("Going to download record %i to %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(
        db="Taxonomy",
        rettype="gb",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=webenv,
        query_key=query_key,
        idtype="acc",
        )
    data = fetch_handle.read()
    myfile= open('rawFile.txt', 'w')
    myfile.writelines(data)
    myfile.close()
    os.system("grep -v ':' rawFile.txt > loading.txt")
    os.system("grep '\.' loading.txt > rawFile.txt")
    os.system("grep -v '\[' rawFile.txt > loading.txt")
    os.system("sed 's/^.*\ //' loading.txt >> GeneList.txt")
    os.rename('GeneList.txt', './GeneLists/{}_GeneList'.format(term))
    os.system('rm ./loading.txt')
    os.system('rm ./rawFile.txt')