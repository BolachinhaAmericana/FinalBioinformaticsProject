#!/usr/bin/env python3
'''
This file has some general functions that will be user throughout the programs
'''
from Bio import Entrez

def entrez_search(database: str, term: str): ######
    '''
    Performs a search query on NCBI

    Arguments:
        Takes the NCBI database type and search term as arguments.
    Vars -
        search: Uses the Entrez API to search our query with history and a retmax of 9999
        search_result: reads eSearch in a variable.
    Returns:
        search_result
    '''
    Entrez.email = 'noWarningPlease@gmail.com'

    if database == 'Nucleotide':
        search = Entrez.esearch(db=database, term=term, usehistory="y",retmax=1,sort='pub date')
    if database == 'Gene':
        search = Entrez.esearch(db=database, term=term, usehistory="y", idtype="acc")
    else:
        search = Entrez.esearch(db=database, term=term, usehistory="y", retmax = 9999)
    search_result = Entrez.read(search)
    return search_result

def entrez_fetch(database: str, web_environment, query_key): ######
    '''
    fetches data using Entrez API

    Arguments:
        takes one NCBI database, the querykey and the web environment from a eSearch
    Vars -
        fetch_handle: Uses Entrez API to fetch information with specified paramethers
        fetch_result: Reads previous value
    Returns:
        fetch_result
    '''
    if database == 'Nucleotide':
        fetch_handle = Entrez.efetch(db=database,
                                     webenv=web_environment,
                                     query_key=query_key,
                                     rettype='fasta',
                                     retmax=1)
    else:
        fetch_handle = Entrez.efetch(db=database, webenv=web_environment, query_key=query_key)
    fetch_result = fetch_handle.read()
    return fetch_result

def print_to_file(output_path, mode, output): ######
    '''
    prints value to file
    Arguments:
        takes a path to output the file and its mode (r, w, a) for example. And the output to print.
    Vars -
        output_file: sets the output file
    Returns:
        None
    '''
    output_file= open(output_path, mode, encoding= 'utf8')
    print(output, file= output_file)
    return None

def file_reader(path: str): ######
    '''
    reads file to var
    Arguments:
        file path
    Vars-
        selected_file: opened file
        file_content: read file
    Return:
        file_content
    '''
    selected_file= open(path, 'r', encoding='utf8')
    file_content= selected_file.read()
    selected_file.close()
    return file_content

def file_reader_to_list(path_to_file: str): ######
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
    with open(path_to_file, 'r', encoding= 'utf8') as file_content_list:
        content_list = [content_list.rstrip() for content_list in file_content_list]
    content_list= list(set(content_list))
    if "" in content_list:
        content_list.remove("")
    return content_list
