#!/bin/python

import sys


# Accepts User Arguments


def getUserArgs():
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

term, taxonTank= getUserArgs()
print(term, taxonTank)

