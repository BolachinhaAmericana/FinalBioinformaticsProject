#!/usr/bin/python
import sys
import taxoniq


def getUserArgs():
    term = sys.argv[1]
    rankTaxonomia = sys.argv[2]
    return term, rankTaxonomia


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


