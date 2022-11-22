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
    ranksList = t.ranked_lineage
    return specieName, commonName, ranksList

if __name__ == '__main__':
    
    term = tID.getUserArguments()
    searchResult = tID.getSearchResult(term) #full eSearch
    taxonId = tID.getTaxonId(searchResult)
    t = loadTaxon(taxonId)
    specieName, commonName, ranksList = getSpecieInfo(t)
    #specieName = getSpecieInfo(t)

    print(specieName)
    print(commonName)
    #print(ranksList)