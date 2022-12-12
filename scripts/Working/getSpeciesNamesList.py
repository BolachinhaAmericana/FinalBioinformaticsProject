#!/bin/python3

import subprocess
import os

from essentials import eSearch, printOutputToFile, readFile, eFetch
import getSpeciesRankName as st


webEnv, queryKey= eSearch('Taxonomy', st.RankName+'[subtree]')["WebEnv"], eSearch('Taxonomy', st.RankName+'[subtree]')["QueryKey"]
#get eFetch.xml
printOutputToFile('./eFetch.xml', 'w', eFetch('Taxonomy', webEnv, queryKey))
# separate all lines to filthering.xml file
filther= readFile('./eFetch.xml').replace('<ScientificName>', '\n')
printOutputToFile('./filthering.xml', 'w', filther)
# Filthers all lines that don't have the word "species"
subprocess.call('./bashScripts/getSpeciesFilthered.sh')
# Separates all line from species name to another line
filther= readFile('./loading.xml').replace('</ScientificName>', '\n')
printOutputToFile('./finalFilther.xml', 'w', filther)
# filthers all lines that are not the species name.
os.system("grep -v 'species' ./finalFilther.xml > ./ScientificNames_list.txt")
os.system('rm ./*.xml')
#subprocess.call('./bashScripts/removeDump.sh')


# Outputs ScientificNames_list.txt