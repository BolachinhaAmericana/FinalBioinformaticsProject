#!/usr/bin/python3
import sys
import taxoniq


def obterArgumentosDoUtilizador():
    term = sys.argv[1]
    rankTaxonomia = sys.argv[2]
    return term, rankTaxonomia


def obterNomeRankOrganismo(term,rankTaxonomia):  
    nomeCientifico = taxoniq.Taxon(scientific_name=term)
    listaRanks = nomeCientifico.ranked_lineage
    relacaoDeRankOrganismo = [(nomeCientifico.rank.name, nomeCientifico) for nomeCientifico in listaRanks]
    for i in relacaoDeRankOrganismo:
        if rankTaxonomia in str(i):
            taxonomia = i[1]
            nomeRankOrganismo = taxonomia.scientific_name
            break
    #print(nomeRankOrganismo)    
    return nomeRankOrganismo    


