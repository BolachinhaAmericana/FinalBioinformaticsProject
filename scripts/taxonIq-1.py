#! /usr/bin/python

import sys

import getTaxonId as tID
import taxoniq

def loadTaxon(taxonId):
    t = taxoniq.Taxon(taxonId)
    return t

def getSpecieInfo(t):
    specieName = t.scientific_name
    commonName = t.common_name
    ranksList = specieName.ranked_lineage
    return specieName, commonName, ranksList

if __name__ == '__main__':
    
    searchResult = tID.getSearchResult("Homo Sapiens") #full eSearch
    taxonId = tID.getTaxonId(searchResult)
    print(taxonId)
    #t = loadTaxon(taxonId)
    #print(t)
    #specieName, commonName, ranksList = getSpecieInfo(t)
    #print(specieName)