#!/usr/bin/env python3

#Made by Valente feat. Silva
from essentials import eSearch, printOutputToFile, readFileToList
from Bio import Entrez
import os


def downloadSpecieGeneList(specie: str, batch_size: int): #4.1
    '''
    Downloads one species gene list
    Arguments:
        specie: scientific name of a specie to get gene list
        batch_size: max= 10000. The number of genes that it will fetch at a time.
    Vars-
        eSearchResults: results of a eSearch on the Gene db, filthering all replaced or discontinued genes.
        count: number of total genes that will be fetched
        webenv: Web Environment
        query_key: Query Key
        fetch_handle: eFetch Query.
        data: read of the fetch_handle
    Returns: None, but will create a LoadingGeneList.txt file

    '''
    eSearchResults= eSearch('Gene', f'{specie} AND alive[prop]')
    count = int(eSearchResults["Count"])
    webenv = eSearchResults["WebEnv"]
    query_key = eSearchResults["QueryKey"]

    for start in range(0, count, batch_size):
        while True:
            try:
                print(f"Downloading record {start+1} to {count} - {specie}")
                # Can't use eFetch function because if has a lot of extra paramethers and might aswell just leave it be.
                fetch_handle = Entrez.efetch(db="Gene",rettype="gb",retmode="text",retstart=start,retmax=batch_size,webenv=webenv,query_key=query_key,idtype="acc",)
                data = fetch_handle.read()
                fetch_handle.close()
                printOutputToFile('loading.dump', 'w', data)
                # for some reason isn't accepting subprocess. So we used os instead ,that works fine.
                #subprocess.call('./bashScripts/getGeneListFilthered.sh')
                
                os.system("grep -v ':' loading.dump > loading1.dump")
                os.system("grep '\.' loading1.dump > loading2.dump")
                os.system("grep -v '\[' loading2.dump > loading3.dump")
                os.system("sed 's/^.*\ //' loading3.dump >> LoadingGeneList.txt")
                os.system('rm ./*.dump')
                break
            except:
                print('Error while downloading. Trying again...')
    return None #Outputs the file  LoadingGeneList.txt

def organizeSpecieGeneList(specie: str, pathToGeneList):
    '''
    Organizes the gene list and reacts incase list is empty.
    Arguments:
        takes the specie name as an argument
    Vars-
        None
    Returns:
        None, but will rename and move a file to GeneLists/
    '''
    try:
        os.rename(pathToGeneList, './GeneLists/{}_GeneList'.format(specie))
    except:
        printOutputToFile('./GeneLists/noGeneSpecies.txt', 'a', specie)
        print(f'The following species has zero genes: {specie}')
    return None
#3- Done- Genes Lists

if __name__ == "__main__":
    speciesList= readFileToList('./ScientificNames_list.txt')
if os.path.exists('./GeneLists'):
    os.system('rm ./GeneLists/*')

else: os.system('mkdir ./GeneLists')
for specie in speciesList:
    downloadSpecieGeneList(specie, 7500)
    organizeSpecieGeneList(specie, './LoadingGeneList.txt')
#3- Done

