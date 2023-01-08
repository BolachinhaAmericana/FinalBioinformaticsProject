#!/usr/bin/env python3

#Made by Silva feat. Marine
import shutil
import sys
import os

def getLoadindDir(geneLists_directory):
    '''
    What for:
        Create a new directory, and copy all the files from the 'geneLists_directory' into this new directory.
        If the 'filteredProximity_GeneLists' already exists, it will delete the directory and try creating it again.
    
    Arguments:
        geneLists_directory: The path to the directory containing the gene lists to be copied.

    Vars:
        directory: string representing the name of the newly created directory.
        filename: string representing the name of the file being copied from the 'geneLists_directory'.

    Returns:
        directory: The name of the newly created directory.

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
    '''
    What for:
        Get the user wanted term, proximity, and similarity values.

    Arguments/Vars:
        term: The user-specified term.
        proximity: The user-specified proximity.
        similarity: The user-specified similarity.

    Returns:
        tuple: A tuple containing the term, proximity, and similarity values.
    '''
    term = sys.argv[1]
    proximity = sys.argv[2]
    similarity = sys.argv[3]
    return term,proximity,similarity

def getsetGoldList(term,directory):
    '''
    What for: 
        Get the gene list for the specified term. 

    Arguments: 
        term: The term for which to get the gene list.
        directory: The directory containing the gene lists.

    Vars:
        listPath: This string represents the path to the file being processed in the for loop. 
        pathToGoldList: This string represents the path to the gene list file.
        geneGoldList: This is a list of strings representing the genes names from the gene list file.
        setGoldList: This is a set of the geneGoldList, representing the gene names from the gene list file.

    Returns:
        A tuple containing the set of genes in the gene list (setGoldList), the directory containing the gene lists (directory), and the path to the gene list file (pathToGoldList).

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
    What for: 
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
    What for:
        Create an empty dictionary with the genes in 'setGoldList' as keys and values initialized to 0.

    Arguments: 
        setGoldList: The set of genes to use as keys in the dictionary.

    Vars:
        goldDict: A dictionary with the genes in 'setGoldList' as keys and values initialized to 0.

    Returns:
        goldDict

    '''
    goldDict = {gene:0 for gene in setGoldList}
    return goldDict


def getGoldDictValue(setGoldList, directory, goldDict, proximity):
    '''
    What for:
        Update the values in the 'goldDict' dictionary for the genes that intersect with the gene lists in the 'directory'.

    Arguments: 
        setGoldList: Set of genes in the gold gene list.
        directory: Directory containing the gene lists.
        goldDict: Dictionary with the genes in 'setGoldList' as keys and values to be updated.
        proximity: Proximity value.

    Vars:
        intresectCount: Integer that keeps track of the number of gene lists that intersect with the gold list.
        listPath: String that represents the file path to a file in the directory.
        geneList: List that represents the genes in a gene list file.
        setGeneList: Represents the gene list read from a file in the directory. 
        gene: Represents a gene in the intersection of the gold list and the gene list read from a file in the directory.

    Returns:
        A tuple containing the updated 'goldDict' dictionary and the count of intersecting gene lists.

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
    '''
    What for:
        Read the names of all the files in a given directory and rewrite them in 'FiltredScientificNames_list.txt'. The names of the files are assumed to be scientific names.

    Arguments: 
        directory: The path to the directory where the scientific names are stored.

    Vars:
        listPath: Represents the file path to a file in the directory.
        scientificName: Represents the scientific name of an organism, extracted from the file.
        filtredScientificNames_list: Represents the file where the scientific names are being written.

    Returns:
        none

    '''
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
    What for: 
        Determine whether the value of a given metric is above a specified threshold.

    Arguments: 
        directory: Path to a directory.
        value: Value of the metric to be compared to the threshold.
        intresectCount: Number of elements in the set being compared.
        similarity: Threshold value to be compared to the metric value.

    Vars:
        none

    Returns:
        True if the value of the metric is above the threshold, False otherwise.

    '''
    similarity = int(similarity)
    if similarity > 100 or similarity < 0:
        shutil.rmtree(directory)
        sys.exit("error: similarity(%) > 100")  
    if value/intresectCount*100 >= similarity:
        return True
    return False    
    

def FiltredGeneNames_list(directory,goldDict,intresectCount,similarity):
    '''
    What for: 
        Write the names of genes that pass a given threshold to a file.

    Arguments: 
        directory: Path to a directory.
        value: Value of the metric to be compared to the threshold.
        goldDict: Dictionary that stores gene names as keys and the number of occurrences as values.
        intresectCount: Number of elements in the set being compared.
        similarity: Threshold value to be compared to the metric value.

    Vars:
        value: Rrepresents the number of occurrences of a gene in the goldDict dictionary.
        filtredGeneNames_list: Represents the file where the gene names that pass the threshold are being written.

    Returns:
        goldDict: Updated dictionary.


    '''
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
    '''
    What for:
        Check if a directory contains only the gold list file. If it does, delete the directory and exit the program. Otherwise, delete the directory.

    Arguments:
        directory: Path to the directory to be checked.
        pathToGoldList: File path to the gold list.

    Vars:
        pathToGoldGene: Represents the name of the gold list file, extracted from pathToGoldList.

    Returns:
        none

    '''
    filenames = [filename for filename in os.listdir(directory)]
    pathToGoldGene = pathToGoldList.replace("filteredProximity_GeneLists/",'')
    if len(filenames) == 1 and pathToGoldGene in filenames:
        shutil.rmtree(directory)
        sys.exit(f"Error: Only input term {pathToGoldGene} have this level of proximity(%) and similarity(%)")
    #Comment the line below to not delete the directory with filtered genes by proximity.    
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