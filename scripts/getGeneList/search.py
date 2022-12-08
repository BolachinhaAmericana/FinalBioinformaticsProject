#!/bin/bash

from Bio import Entrez
import sys

Entrez.email = 'duarte.tz.valente@gmail.com'


def getUserArgs():

    db = sys.argv[1]
    term = sys.argv[2]
    return db, term

def eSearch(db, term):

    eSearch = Entrez.esearch(db=db, term=term, usehistory="y", retmax = 10)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

def getWebEnv(eSearch):

    webEnv = eSearch["WebEnv"]
    return webEnv

def getQueryKey(eSearch):

    queryKey = eSearch["QueryKey"]
    return queryKey

def eFetch(queryKey, webEnv):

    fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey, rettype='txt', lvl = 20)
    eFetchResult = fetchHandle.read()
    return eFetchResult

if __name__ == '__main__':
    #db, term = obterArgumentosDoUtilizador()
    db = "gene"
    #term = "Passer domesticus"
    #term = "Canis lupus"
    term = 'Panthera leo'
    eSearchResult = eSearch(db, term)
    queryKey, webEnv = getQueryKey(eSearchResult), getWebEnv(eSearchResult)
    eFetchResult = eFetch(queryKey, webEnv)
    with open('./eSearch.xml', 'w') as sys.stdout:
        print(eSearchResult)
    with open('./eFetch.xml', 'w') as sys.stdout:
        print(eFetchResult)
