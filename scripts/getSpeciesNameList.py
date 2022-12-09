#!/bin/python

from Bio import Entrez
import sys
import taxoniq
import os

Entrez.email = 'NoWarningPls@gmail.com'

#term = 'Homo sapiens'
#taxonRank = 'order'

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

term, taxonRank= getUserArgs()

def eSearch(db, term):
    '''
    Performs a search query on NCBI 

    Arguments: 
        Takes the NCBI database type and search term as arguments.
    Vars -
        eSearch: Uses the Entrez API to search our query with history and a retmax of 9999
        eSearchResult: reads eSearch in a variable.
    Returns:
        eSearchResult
    '''
    eSearch = Entrez.esearch(db=db, term=term, usehistory="y", retmax = 9999)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

eSearchResult= eSearch('Taxonomy', term)

def getTaxonId(eSearchResult):
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

taxonId= getTaxonId(eSearchResult)

def getNameRank(taxonId, taxonRank):
    '''
    Will get the specific name of the rank according to species and overall rank. Example: Human, order = Primates

    Arguments:
        Takes target taxon id and input taxon rank as arguments.
    Vars -
        t: loads taxon on var using Taxoniq API
        rankList: List of lineage of loaded Taxon
        relation: Compares the lineage with taxon ranks
        tax: taxon id of the related rank
        taxNameRank: conversion of taxon id to scientific name using Taxoniq API
    Returns:
        taxNameRank
    '''
    t = taxoniq.Taxon(taxonId)
    rankList = t.ranked_lineage
    relation = [(t.rank.name, t) for t in rankList]
    for i in relation:
        if taxonRank in str(i):
            tax = i[1]
            taxNameRank = tax.scientific_name
            break
    return taxNameRank

taxNamerank= getNameRank(taxonId, taxonRank)

eSearchResult= eSearch('Taxonomy', taxNamerank+'[subtree]')

def getWebEnv(eSearch):
    '''
    gets eSearch Web Environment

    Arguments: 
        takes eSearch Result as argument
    Vars - 
        webEnv: contain web environment
    Returns:
        webEnv
    '''
    webEnv = eSearch["WebEnv"]
    return webEnv
def getQueryKey(eSearch):
    '''
    gets eSearch Query Key

    Arguments: 
        takes eSearch Result as argument
    Vars - 
        queryKey: contain Query Key
    Returns:
        queryKey
    '''
    queryKey = eSearch["QueryKey"]
    return queryKey
webEnv= getWebEnv(eSearchResult)
queryKey= getQueryKey(eSearchResult)

def eFetch(db, queryKey, webEnv):
    '''
    fetches data using Entrez API

    Arguments:
        takes one NCBI database, the querykey and the web environment from a eSearch
    Vars -
        fetchHandle: Uses Entrez API to fetch information with specified paramethers
        eFetchResult: Reads previous value
    Returns:
        eFetchResult
    '''
    if db== 'Taxonomy':
        #fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey, rettype='uilist', retmode= 'text')
        fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey)
        eFetchResult = fetchHandle.read()
    return eFetchResult

def downloadEFetch():
    '''
    downloads eFetch search as .xml

    Arguments:
        None
    Vars- 
        record: executes eFetch function instance to download
    Return:
        None, but wull create a eFetch.xml file
    '''
    record= eFetch('Taxonomy', queryKey, webEnv)
    with open('./eFetch.xml', 'w') as sys.stdout:
            print(record)
    return None

downloadEFetch()

def getSpeciesFilthered():
    '''
    Filthers all instances that are not a specie in eFetch.xml file

    Arguments:
        None
    Vars -
        f: opens file
        file: reads file
        first: filthers file
    Return:
        None, but will create a loading.xml file
    '''
    f = open('eFetch.xml', 'r')
    file= f.read()
    f.close()
    first= file.replace('<ScientificName>', '\n')
    with open('test.xml', 'w') as sys.stdout:
        print(first)
    os.system("grep 'species' test.xml > list.xml")
    os.system("grep -v 'subspecies' list.xml > loading.xml")
    os.system('rm list.xml')
    os.system('rm test.xml')
    return None

getSpeciesFilthered()

def getSpeciesNameList():
    '''
    Filthers loading.xml to create a file with all scientific names of the species in study.

    Arguments:
        None
    Vars- 
        f: opens file
        file: reads file
        first: filthers file
    Returns: 
        None, but will create a ScientificNames_list.txt file    
    '''

    f = open('loading.xml', 'r')
    file= f.read()

    first= file.replace('</ScientificName>', '\n')

    with open('test1.xml', 'w') as sys.stdout:
        print(first)

    os.system("grep -v 'species' test1.xml > ScientificNames_list.txt")
    os.system('rm loading.xml')
    os.system('rm test1.xml')
    os.system('rm eFetch.xml')
    return None

getSpeciesNameList()