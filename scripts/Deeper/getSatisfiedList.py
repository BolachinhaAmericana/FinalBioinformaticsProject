#!/bin/python3
import shutil
import sys
import os

def getLoadindDir(geneLists_directory):
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
    proximity = int(proximity)
    if proximity > 100:
        shutil.rmtree(directory)
        sys.exit("error: proximity(%) > 100")
    if int(len(setGoldList.intersection(setGeneList))/len(setGoldList)*100) >= proximity:
        return True
    return False    


def getGoldEmptyDict(setGoldList):
    goldDict = {gene:0 for gene in setGoldList}
    return goldDict


def getGoldDictValue(setGoldList, directory, goldDict,proximity):
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
            os.remove("ScientificNames_list.txt")
            filtredScientificNames_list = open('FiltredScientificNames_list.txt', "x")
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
    similarity = int(similarity)
    if similarity > 100:
        shutil.rmtree(directory)
        sys.exit("error: similarity(%) > 100")  
    if value/intresectCount*100 >= similarity:
        return True
    return False    


def getsimilarityGoldDict(directory,goldDict,intresectCount,similarity):
    for gene in list(goldDict.keys()):
        value = goldDict[gene]
        if(not(haveSimilarity(directory,value,intresectCount,similarity))):
            del goldDict[gene]
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
    similarityGoldDict = getsimilarityGoldDict(directory,goldDict,intresectCount,similarity)
    #print(similarityGoldDict)
    verifyDir(directory,pathToGoldList)
    
       

