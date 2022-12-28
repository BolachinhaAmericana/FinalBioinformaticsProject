#!/usr/bin/env python3

#Made by Valente feat. Silva
#V.Silva

import sys
from Bio import Entrez

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

    if db == 'Nucleotide':
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y",retmax=1,sort='pub date')
    if db == 'Gene':
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y", idtype="acc")
    else:
        eSearch = Entrez.esearch(db=db, term=term, usehistory="y", retmax = 9999)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

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
    if db == 'Nucleotide':
        fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey, rettype='fasta',retmax=1)
    else:
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
