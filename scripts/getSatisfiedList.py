#!/usr/bin/env python3

#Made by Silva feat. Marine
import shutil
import sys
import os

def getLoadindDir(geneLists_directory):
    '''
    What for 
		Tells where to create the directory 
        Creates an temporary folder with a gene list

    Arguments: 
        

    Vars:
         

    Returns:
        

    '''
    while True:
        try:
            os.mkdir('filteredProximity_GeneLists')
            directory = 'filteredProximity_GeneLists'
            for filename in os.listdir(geneLists_directory):
                shutil.copy2(os.path.join(geneLists_directory,filename), directory)
            return directory
        except FileExistsError:
            shutil.rmtree('filteredProximity_GeneLists')    

def getUserArguments(term,proximity,similarity):
    term = sys.argv[1]
    proximity = sys.argv[2]
    similarity = sys.argv[3]
    return term,proximity,similarity

def getsetGoldList(term,directory):
    '''
    What for 
		Selects a directory and a term 

    Arguments: 
        

    Vars:
         

    Returns:
        Directory with the gold list and the gold list

    '''
    try:
        for filename in os.listdir(directory):
            listPath = os.path.join(directory, filename)
            pathToGoldList = f"{directory}/{term}_GeneList"
            if listPath == pathToGoldList:
                with open(pathToGoldList, 'r') as goldList:
                    geneGoldList = [nameList.rstrip() for nameList in goldList]
                setGoldList = set(geneGoldList)  
        return setGoldList, directory ,pathToGoldList
    except UnboundLocalError:
        print("Bad Specie name")


def intersect(directory,setGeneList,setGoldList,proximity):
    '''
    What for 
		Calculates the proximity percentage of the Gold gene related to the genes of other organisms of the list

    Arguments: 
        directory: wanted directory
        setGeneList: List of genes
        setGoldList: Animal of input
        proximity: proximity percentage

    Vars:
        none

    Returns:
        logical value if the proximity percentage calculated is equal or superior of the proximity percentage inputed

    '''
    proximity = int(proximity)
    if proximity > 100 or proximity < 0:
        shutil.rmtree(directory)
        sys.exit("error: proximity(%) > 100")
    if int(len(setGoldList.intersection(setGeneList))/len(setGoldList)*100) >= proximity:
        return True
    return False    


def getGoldEmptyDict(setGoldList):
    '''
    What for 
		 

    Arguments: 
        

    Vars:
        none

    Returns:
        goldDict

    '''
    goldDict = {gene:0 for gene in setGoldList}
    return goldDict


def getGoldDictValue(setGoldList, directory, goldDict,proximity):
    '''
    What for 
		cria novo conjunto e genes para cada gene 
        se nao de remove 
        se der adiciona +1 gene em todos os animas que for encontrado 

    Arguments: 
        

    Vars:
         

    Returns:
        True or false value 

    '''
    intresectCount = 0
    for filename in os.listdir(directory):
        listPath = os.path.join(directory, filename)
        with open(listPath, 'r') as geneList:
            geneList = [nameList.rstrip() for nameList in geneList]
        setGeneList = set(geneList)
        if not(intersect(directory,setGeneList, setGoldList,proximity)):
            os.remove(listPath)
            continue
        intresectCount += 1    
        setGeneList = setGoldList.intersection(setGeneList)
        for gene in setGeneList:
            goldDict[gene] += 1
    return  goldDict, intresectCount          

def getFilteredScientificName_list(directory):
    while True:
        try:
            #os.remove("ScientificNames_list.txt")
            for filename in os.listdir(directory):
                listPath = os.path.join(directory, filename)
                listPath = listPath.replace("filteredProximity_GeneLists/",'')
                scientificName = listPath.replace("_GeneList",'')   
                filtredScientificNames_list = open("FiltredScientificNames_list.txt", "a")
                filtredScientificNames_list.writelines("%s\n" %scientificName)
                filtredScientificNames_list.close()
            break
        except FileNotFoundError:
            open('ScientificNames_list.txt', "x")
        except FileExistsError:
            os.remove("FiltredScientificNames_list.txt")           

def haveSimilarity(directory,value,intresectCount,similarity):
    '''
    What for 
		Tests the similarity 

    Arguments: 
        

    Vars:
         

    Returns:
        

    '''
    similarity = int(similarity)
    if similarity > 100 or similarity < 0:
        shutil.rmtree(directory)
        sys.exit("error: similarity(%) > 100")  
    if value/intresectCount*100 >= similarity:
        return True
    return False    
    

def FiltredGeneNames_list(directory,goldDict,intresectCount,similarity):
    while True:
        try:
            os.remove("FiltredGeneNames_list.txt")
            for gene in list(goldDict.keys()):
                value = goldDict[gene]
                if(not(haveSimilarity(directory,value,intresectCount,similarity))):
                    del goldDict[gene]
                    continue
                filtredGeneNames_list = open("FiltredGeneNames_list.txt", "a")
                filtredGeneNames_list.writelines("%s\n" %gene)
                filtredGeneNames_list.close()
            break
        except FileNotFoundError:
            open('FiltredGeneNames_list.txt', "x")
    return goldDict       



def verifyDir(directory,pathToGoldList):
    filenames = [filename for filename in os.listdir(directory)]
    pathToGoldGene = pathToGoldList.replace("filteredProximity_GeneLists/",'')
    if len(filenames) == 1 and pathToGoldGene in filenames:
        shutil.rmtree(directory)
        sys.exit(f"Error: Only input term {pathToGoldGene} have this level of proximity(%) and similarity(%)")
    #Comentar linha a baixo para não apagar dir com genes filtrados por proximidade    
    shutil.rmtree(directory)

if __name__ == "__main__":
    directory = getLoadindDir(geneLists_directory='GeneLists')
    term,proximity,similarity = getUserArguments(term='',proximity='',similarity='')
    setGoldList, directory, pathToGoldList = getsetGoldList(term,directory)
    goldDict = getGoldEmptyDict(setGoldList)
    goldDict,intresectCount = getGoldDictValue(setGoldList, directory, goldDict,proximity)
    getFilteredScientificName_list(directory)
    filtredDict = FiltredGeneNames_list(directory,goldDict,intresectCount,similarity)
    #print(filtredDict)
    verifyDir(directory,pathToGoldList)