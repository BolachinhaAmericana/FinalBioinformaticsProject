#!/usr/bin/python3
import sys
import taxoniq


def getArgumentFromUser():
    term = sys.argv[1]
    rankTaxonomy = sys.argv[2]
    return term, rankTaxonomy 


def getNameRankOrganism(term,rankTaxonomy):
    nameScientific  = taxoniq.Taxon(scientific_name=term)
    listRanks = nameScientific.ranked_lineage
    relactionRankOrganism = [(nameScientific.rank.name, nameScientific) for nameScientific in listRanks]
    for i in relactionRankOrganism:
        if rankTaxonomy in str(i):
            taxonomy = i[1]
            nameRankOrganism = taxonomy.scientific_name
            break
    #print(nameRankOrganism)
    return nameRankOrganism