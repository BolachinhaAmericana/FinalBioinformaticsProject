import os

def getFiltheredSpeciesNameList(pathToFiltheredLists):
    if os.path.exists('FiltredNamesList.txt'):
        os.remove('FiltredNamesList.txt')

    for filename in os.listdir(pathToFiltheredLists):
        listPath= os.path.join(pathToFiltheredLists, filename)
        listPath= listPath.replace(f"{pathToFiltheredLists}/",'')
        nameList= listPath.replace("_GeneList", "")
        filtheredNamesList = open("FiltredNamesList.txt", "a")
        filtheredNamesList.writelines(nameList+'\n')
        filtheredNamesList.close()

getFiltheredSpeciesNameList('FiltheredLists')