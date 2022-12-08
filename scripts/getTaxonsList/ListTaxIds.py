#!/bin/bash

from Bio import Entrez
import sys
import os
import taxoniq

def getUserArgs():
    '''
    Gets user arguments: 
        term = search term (this will always be a specie with only 1st letter as Capital)
        rankTaxonomy = desired taxonomy rank (options: )
    
    '''
    term = sys.argv[1]
    rankTaxonomy = sys.argv[2]
    return term, rankTaxonomy
    
def getSearchResult(nomeRankOrganismo):
    eSearch = Entrez.esearch(db="Taxonomy", term=nomeRankOrganismo, usehistory='y', retmax=20)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

def getTaxonId(eSearchResult):
    taxonId = str(eSearchResult["IdList"])
    taxonId = taxonId.replace("[\'",'')
    taxonId = taxonId.replace("\']",'')
    return taxonId

def getNameRank(taxonId, taxonRank):
        t = taxoniq.Taxon(taxonId)
        rankList = t.ranked_lineage
        relation = [(t.rank.name, t) for t in rankList]
        for i in relation:
            if taxonRank in str(i):
                tax = i[1]
                taxNameRank = tax.scientific_name
                break
        return taxNameRank

def getRelatedTaxIdList():
    term, rankTaxonomia = getUserArgs()
    searchResults = getSearchResult(term)
    taxonId = getTaxonId(searchResults)
    nameRank = getNameRank(taxonId, rankTaxonomia)
    rankSearchResults = getSearchResult(nameRank)
    rankTaxonId = getTaxonId(rankSearchResults)

    os.system(f"taxonkit list --ids {rankTaxonId} -nr --indent '' > loading")
    os.system(f"sed -e 's/sub*.//' loading | grep -e 'species' | sed -e 's/.*species]//' | sed -e 's/ //' > namesByRank.txt")
    os.system(f"rm -r loading")
    return None

if __name__ == '__main__':
    term, rankTaxonomy = getUserArgs()
    searchResults = getSearchResult(term)
    taxonId = getTaxonId(searchResults)
    nameRank = getNameRank(taxonId, rankTaxonomy)
    getRelatedTaxIdList()




