import essentials as e
import 3-

import os

def organizeSpecieGeneList(specie: str, pathTogeneList):
    '''
    Organizes the gene list and reacts incase list is empty.
    Arguments:
        takes the specie name as an argument
    Vars-
        None
    Returns:
        None
    '''
    try:
        os.rename(pathTogeneList, './GeneLists/{}_GeneList'.format(specie))
    except:
        e.printOutputToFile('./GeneLists/noGeneSpecies.txt', 'a', specie+'\n')
        print(f'The following species has zero genes: {specie}')
    return None

speciesList= e.readFileToList('./ScientificNames_list.txt')

if os.path.exists('./GeneLists'):
    os.system('rm ./GeneLists/*')
else: os.system('mkdir ./GeneLists')

for specie in speciesList:
    rd.downloadSpecieGeneList(specie, 7500)
    organizeSpecieGeneList(specie, './LoadingGeneList.txt')
