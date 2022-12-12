import subprocess
import os

import essentials as e
import 1-


webEnv, queryKey= e.eSearch('Taxonomy', st.RankName+'[subtree]')["WebEnv"], e.eSearch('Taxonomy', st.RankName+'[subtree]')["QueryKey"]
#get eFetch.xml
e.printOutputToFile('./eFetch.xml', 'w', e.eFetch('Taxonomy', webEnv, queryKey))
# separate all lines to filthering.xml file
filther= e.readFile('./eFetch.xml').replace('<ScientificName>', '\n')
e.printOutputToFile('./filthering.xml', 'w', filther)
# Filthers all lines that don't have the word "species"
subprocess.call('./bashScripts/getSpeciesFilthered.sh')
# Separates all line from species name to another line
filther= e.readFile('./loading.xml').replace('</ScientificName>', '\n')
e.printOutputToFile('./finalFilther.xml', 'w', filther)
# filthers all lines that are not the species name.
os.system("grep -v 'species' ./finalFilther.xml > ./ScientificNames_list.txt")
os.system('rm ./*.xml')
#subprocess.call('./bashScripts/removeDump.sh')