#!/usr/bin/python


# Get TaxonId by species Name.


import sys
from Bio import Entrez
import os
import subprocess

def getUserArguments():
    term = sys.argv[1]
    return term

def getSearchResult(term):
    eSearch = Entrez.esearch(db="Taxonomy", term=term, usehistory='y', retmax=20)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

def getTaxonId(result):
    taxonId = result["IdList"]
    os.system(f"echo {taxonId} > 1.txt")
    subprocess.call("./clearTaxonId.sh")
    with open('3.txt', 'r') as file:
        clearedTaxonId = file.read().replace('\n', '')
    os.system("rm 3.txt")
    return clearedTaxonId

if __name__ == '__main__':

    term = getUserArguments()
    #term = "Homo Sapiens"
    searchResult = getSearchResult(term)
    taxonId = getTaxonId(searchResult)
    print(taxonId)

