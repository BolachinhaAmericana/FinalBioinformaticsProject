#!/usr/bin/python

import os

import getTaxonId
import getRankName

def run():
    term, rankTaxonomia = getRankName.getUserArgs()
    searchResults = getTaxonId.getSearchResult(term)
    taxonId = getTaxonId.getTaxonId(searchResults)
    nameRank = getRankName.getNameRank(taxonId, rankTaxonomia)
    rankSearchResults = getTaxonId.getSearchResult(nameRank)
    rankTaxonId = getTaxonId.getTaxonId(rankSearchResults)

    os.system(f"taxonkit list --ids {rankTaxonId} -nr --indent '' > output.txt")
    os.system(f"grep -r 'species' ./output.txt > 1.txt")
    os.system(f"grep -v 'subspecies' ./1.txt > output.txt")
    os.system(f"sed 's/\s.*$//' output.txt  > ListTaxonIds.txt")
    os.system(f"rm 1.txt")
    os.system(f"rm output.txt")
    return None


if __name__ == '__main__':
    run()