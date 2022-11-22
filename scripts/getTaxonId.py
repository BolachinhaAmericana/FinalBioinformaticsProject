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
    #print(result)
    os.system(f"echo {taxonId} > 1.txt")
    clearedTaxonId = 'sim'
    subprocess.call("./clearTaxonId.sh")
    return

if __name__ == '__main__':

    term = getUserArguments()
    searchResult = getSearchResult(term)
    getTaxonId(searchResult)
