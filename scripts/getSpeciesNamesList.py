#!/usr/bin/env python3
'''
This file will output a file that contains the list of every species related
Outputs ScientificNames_list.txt
'''

import sys
import os
from essentials import entrez_search, entrez_fetch, print_to_file, file_reader

def get_user_arguments(rank_name):
    '''
    What for:
        Retrieves takes arguments passed by the user.

    Arguments:
        None

    Vars:
        rank_name: The rank name input by the user.

    Returns:
        rank_name
    '''
    try:
        rank_name = sys.argv[1]
        return rank_name
    except IndexError:
        sys.exit('No input was give. This program takes a species taxonomy level (ex: primates).')

def scientific_names_list_downloader(rank_name: str):
    '''
    What for:
        Downloads ScientificNames_list.txt, with all species from said taxonomy rank name.

    Arguments:
        rank_name: Name of the taxonomy rank according to said specie.

    Vars:
        web_environment and query_key: Web Environment and Query Key
        filther: filthering instance from efetch query.

    Returns:
        None
    '''
    web_environment= entrez_search('Taxonomy', rank_name+'[subtree]')["WebEnv"]
    query_key= entrez_search('Taxonomy', rank_name+'[subtree]')["QueryKey"]

    print_to_file('./eFetch.dump', 'w', entrez_fetch('Taxonomy', web_environment, query_key))
    filther= file_reader('./eFetch.dump').replace('<ScientificName>', '\n')
    print_to_file('./filthering.dump', 'w', filther)
    os.system("grep 'species' filthering.dump > list.dump")
    os.system("grep -v 'subspecies' list.dump > loading.dump")
    filther= file_reader('./loading.dump').replace('</ScientificName>', '\n')
    print_to_file('./finalFilther.dump', 'w', filther)

    os.system("grep -v 'species' ./finalFilther.dump > ./ScientificNames_list.txt")
    os.system('rm ./*.dump')
    return None



if __name__ == "__main__":
    RANK_NAME = get_user_arguments(rank_name='')
    scientific_names_list_downloader(RANK_NAME)
