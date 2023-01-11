#!/usr/bin/env python3

#Made by Valente
from essentials import eSearch, printOutputToFile, readFile, eFetch
import sys
import os


def getUserArgs(RankName):

    RankName = sys.argv[1]
    return RankName


def downloadScientificNamesList(RankName: str):
    '''
    Downloads ScientificNames_list.txt, with all species from said taxonomy rank name
    Args:
        RankName: Name of the taxonomy rank according to said specie. Example: Human order = Primates
    Vars-
        webEnv and queryKey: Web Environment and Query Key
        filther: filthering instance from efetch query.
    Returns: 
        None, but will create a ScientificNames_list.txt file with the scientific name of all species from RankName Argument.
    '''
    webEnv, queryKey= eSearch('Taxonomy', RankName+'[subtree]')["WebEnv"], eSearch('Taxonomy', RankName+'[subtree]')["QueryKey"]
    printOutputToFile('./eFetch.dump', 'w', eFetch('Taxonomy', webEnv, queryKey)) #get eFetch
    filther= readFile('./eFetch.dump').replace('<ScientificName>', '\n')
    printOutputToFile('./filthering.dump', 'w', filther)
    os.system("grep 'species' filthering.dump > list.dump")
    os.system("grep -v 'subspecies' list.dump > loading.dump")
    filther= readFile('./loading.dump').replace('</ScientificName>', '\n')
    printOutputToFile('./finalFilther.dump', 'w', filther)
    os.system("grep -v 'species' ./finalFilther.dump > ./ScientificNames_list.txt")
    os.system('rm ./*.dump')
    return None



if __name__ == "__main__":
    RankName = getUserArgs(RankName='') 
    downloadScientificNamesList(RankName)

# Outputs ScientificNames_list.txt
