#!/usr/bin/python

# Get TaxonId.

from Bio import Entrez
import getRankName


def getSearchResult(nomeRankOrganismo):
    eSearch = Entrez.esearch(db="Taxonomy", term=nomeRankOrganismo, usehistory='y', retmax=20)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

def getTaxonId(eSearchResult):
    taxonId = str(eSearchResult["IdList"])
    taxonId = taxonId.replace("[\'",'')
    taxonId = taxonId.replace("\']",'')
    return taxonId

if __name__ == '__main__':
    term, rankTaxonomia = getRankName.getUserArgs()
    searchResults = getSearchResult(term)
    taxonId = getTaxonId(searchResults)
    nameRank = getRankName.getNameRank(taxonId, rankTaxonomia)

    rankSearchResults = getSearchResult(nameRank)
    rankTaxonId = getTaxonId(rankSearchResults)

    print('Specie Taxon Id')
    print(taxonId)
    print('Taxon Rank Name:')
    print(nameRank)
    print('Taxon Rank Id:')
    print(rankTaxonId)


