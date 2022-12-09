from Bio import Entrez
import os
import sys

Entrez.email = "guilherme.vcc.silva@gmail.com" 

def eSearch(db, term):
    '''
    Performs a search query on NCBI 

    Arguments: 
        Takes the NCBI database type and search term as arguments.
    Vars -
        eSearch: Uses the Entrez API to search our query with history and a retmax of 9999
        eSearchResult: reads eSearch in a variable.
    Returns:
        eSearchResult
    '''
    if db=='Gene':
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y", idtype="acc")
    else: 
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y", retmax = 9999)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

term= 'Panthera leo' #One specie at a Time

search_results= eSearch("Gene", term)

acc_list = search_results["IdList"]
count = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]
batch_size = 9999

os.system("rm rawFile.txt")

for start in range(0, count, batch_size):
    end = count
    print("Going to download record %i to %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(
        db="Gene",
        rettype="gb",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=webenv,
        query_key=query_key,
        idtype="acc",
        )
    data = fetch_handle.read()

    myfile= open('rawFile.txt', 'a')
    myfile.writelines(data)
    myfile.close()

    fetch_handle.close()

os.system("grep -v ':' rawFile.txt > loading.txt")
os.system("grep '\.' loading.txt > rawFile.txt")
os.system("grep -v '\[' rawFile.txt > loading.txt")
os.system("sed 's/^.*\ //' loading.txt > GeneList.txt")
os.system("rm loading.txt")
os.system("rm rawFile.txt")
os.rename('GeneList.txt', '{}_GeneList'.format(term))