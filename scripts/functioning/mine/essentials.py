#!/bin/python3

import sys
from Bio import Entrez

def getUserArgs(): #1.1
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

def eSearch(db: str, term: str): #1.2
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
    Entrez.email = 'noWarningPlease@gmail.com'
    if db== 'Gene':
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y", idtype="acc")
    else:
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y", retmax = 9999)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

#term, taxonRank= getUserArgs()
#eSearchResult= eSearch('Taxonomy', term)
#webEnv = eSearchResult["WebEnv"]      # 1.3
#queryKey = eSearchResult["QueryKey"]  # 1.3

def eFetch(db: str, webEnv, queryKey): #1.4
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
    fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey)
    eFetchResult = fetchHandle.read()
    return eFetchResult

def printOutputToFile(path: str, type: str, output): # 1.5
    '''
    temporarily redirects stdout to a file.
    Arguments:
        takes a path to output the file and its type (r, w, a) for example. And the output to print.
    Vars -
        stdout: saves unchanged standard output value
    Returns:
        None
    '''
    stdout= sys.stdout
    sys.stdout= open(path, type)
    print(output)
    sys.stdout= stdout
    return None

def readFile(path: str): #1.6
    '''
    reads file to var
    Arguments:
        file path
    Vars- 
        f: opened file
        file: read file
    Return:
        file
    '''
    f=open(path)
    file= f.read()
    f.close()
    return file

def readFileToList(path_to_file: str): #1.7
    '''
    reads a file line by line. Adds lines to a list and removes repeated and empty items.
    Args:
        takes the path to the file as an argument
    Vars - 
        fileList: unchanged file content
        List: transformation of the file content into a python list.
    Returns:
        cleared and transformed list.
    '''
    with open(path_to_file, 'r') as fileList:
        List = [List.rstrip() for List in fileList]
    List= list(set(List)) 
    if "" in List:
        List.remove("")
    return List

