#!/usr/bin/python


# Get TaxonId by species Name.


import sys
from Bio import Entrez
import os

def getUserArguments():
    term = sys.argv[1]
    return term

def getSearchResult(term):
    eSearch = Entrez.esearch(db="Taxonomy", term=term, usehistory='y', retmax=20)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

def getTaxonId(result):
    taxonId = result["IdList"]
    os.system(f"echo {taxonId} > temp1.txt")
    os.system("sed -i 's/\[//g' temp1.txt >> temp1.txt")
    os.system("sed -i 's/\]//g' temp1.txt")
    clearedTaxonId = os.system("cat temp1.txt")
    os.system("rm temp1.txt")

    return clearedTaxonId


if __name__ == '__main__':

    term = getUserArguments()
    searchResult = getSearchResult(term)
    taxonId = getTaxonId(searchResult)
    print(taxonId)

