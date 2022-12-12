import taxoniq

# Accepts a Taxon Id and a input rank and returns de name of the rank acording to the ID

def getRankName(taxonId, taxonRank):
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

'''
Example output:
Primates
'''