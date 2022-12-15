#!/bin/python3

import os
from Bio import Entrez
#import subprocess

from essentials import eSearch, printOutputToFile 

def downloadSpecieGeneList(specie: str, batch_size: int): #4.1
    '''
    Downloads 
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
    return None #Outputs the file GeneList.txt

#downloadSpecieGeneList('Passer domesticus', 7500)