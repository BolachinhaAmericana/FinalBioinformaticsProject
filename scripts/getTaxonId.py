#!/usr/bin/python3

# Get TaxonId.

from Bio import Entrez
import getRankName as nomeRank

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
    term, rankTaxonomia = nomeRank.obterArgumentosDoUtilizador()
    nomeRankOrganismo = nomeRank.obterNomeRankOrganismo(term,rankTaxonomia)             
    searchResult = getSearchResult(nomeRankOrganismo)
    taxonId = getTaxonId(searchResult)
    print(taxonId)

