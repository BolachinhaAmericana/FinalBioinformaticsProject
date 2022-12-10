from Bio import Entrez
import os
import sys

Entrez.email = "noWarningPls@gmail.com" 

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
        eSearch = Entrez.esearch(db=db, term=f'{term} AND alive[prop]', usehistory="y", idtype="acc")
    else: 
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y", retmax = 9999)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

if os.path.exists('./GeneLists') == False:
    os.system("mkdir ./GeneLists")

with open('ScientificNames_list.txt', 'r') as List:
    nameList = [nameList.rstrip() for nameList in List]

nameList= list(set(nameList)) 
if "" in nameList:
    nameList.remove("")



for name in nameList:

    term= (f'{name}') 

    search_results= eSearch("Gene", term)

    acc_list = search_results["IdList"]
    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    batch_size = 5000

    if os.path.isfile('rawFile.txt')== True:
        os.system('rm rawFile.txt')

    
    for start in range(0, count, batch_size):
        end = count
        print(f"Going to download record {start+1} to {end} - {term}")
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

        myfile= open('rawFile.txt', 'w')
        myfile.writelines(data)
        myfile.close()
        os.system("grep -v ':' rawFile.txt > loading.txt")
        os.system("grep '\.' loading.txt > rawFile.txt")
        os.system("grep -v '\[' rawFile.txt > loading.txt")
        os.system("sed 's/^.*\ //' loading.txt >> GeneList.txt")
        

        fetch_handle.close()

    
    os.rename('GeneList.txt', './GeneLists/{}_GeneList'.format(term))

