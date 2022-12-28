#!/usr/bin/env python3

#Made by Valente feat. Silva
from essentials import eSearch
import taxoniq
import sys


def getUserArgs(term,taxonRank):
    '''
    Takes 2 arguments as std.input
    
    Arguments:
         None
    Vars -
        term: Search term. Should be a species scientific name following the example: "Canis lupus". Note that the 1st letter is the only capital.
        taxonRank: Taxonomic Rank relation. If you dont know what to put here just enter an item of the following list: 
            ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    Returns:
        term and taxonRank
    '''
    term = sys.argv[1]
    taxonRank= sys.argv[2] 
    return term, taxonRank

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

if __name__ == "__main__":
    term, taxonRank= getUserArgs(term='', taxonRank='')
    RankName= getRankName(getTaxonId(eSearch('Taxonomy', term)), taxonRank)
    print(RankName)
