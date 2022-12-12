import essentials as e
import taxoniq

def getTaxonId(eSearchResult): # 2.1
    '''
    Will get the Taxon Id of our input search term

    Arguments: 
        Takes eSearch result as argument.
    Vars - 
        taxonId: filthering string
    Returns:
        Filthered Taxon Id
    '''
    taxonId = str(eSearchResult["IdList"])
    taxonId = taxonId.replace("[\'",'')
    taxonId = taxonId.replace("\']",'')
    return taxonId

def getRankName(taxonId, taxonRank): #2.2
    '''
    Will get the specific name of the rank according to species and overall rank. Example: Human, order = Primates

    Arguments:
        Takes target taxon id and input taxonomy rank as arguments.
    Vars -
        t: loads taxon on var using Taxoniq API
        rankList: List of lineage of loaded Taxon
        relation: Compares the lineage with taxon ranks
        tax: taxon id of the related rank
        taxRankName: conversion of taxon id to scientific name using Taxoniq API
    Returns:
        taxRankName
    '''
    t = taxoniq.Taxon(taxonId)
    rankList = t.ranked_lineage
    relation = [(t.rank.name, t) for t in rankList]
    for i in relation:
        if taxonRank in str(i):
            tax = i[1]
            taxRankName = tax.scientific_name
            break
    return taxRankName

term, taxonRank= e.getUserArgs()

RankName= getRankName(getTaxonId(e.eSearch('Taxonomy', term)), taxonRank)
#print(RankName) # Module Working