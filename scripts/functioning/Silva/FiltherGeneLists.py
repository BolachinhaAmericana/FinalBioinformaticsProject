import sys
import os
import shutil
import essentials as e

def getUserArguments():
    term = sys.argv[1]
    taxonRank = sys.argv[2]
    proximity = sys.argv[3]
    similarity = sys.argv[4]
    return term,taxonRank, proximity,similarity

def getLoadindDir(geneLists_directory):
    unfiltheredDir= geneLists_directory
    while True:
        try:
            os.mkdir('FiltheredGeneLists')
            directory = 'FiltheredGeneLists'
            for filename in os.listdir(geneLists_directory):
                shutil.copy2(os.path.join(geneLists_directory,filename), directory)
            return directory, unfiltheredDir
        except FileExistsError:
            shutil.rmtree('FiltheredGeneLists')

def setGoldenList(term):
    try:
        for file in os.listdir('GeneLists'):
                geneListPath = os.path.join('GeneLists', f"{file}")
                goldListPath= os.path.join('GeneLists', f"{term}_GeneList")
                if geneListPath==goldListPath:
                    with open(goldListPath, 'r') as goldenList:
                        rawGoldenGenes= [nameList.rstrip() for nameList in goldenList]
                    rawGoldenList = set(rawGoldenGenes)
        return rawGoldenList, goldListPath
    except UnboundLocalError:
        print("Bad Specie name")

def speciesValidation(geneList, goldenList, percentage: int): #Filthering species.
    if percentage > 100 or percentage < 0:
            sys.exit("error: proximity(%) > 100")
    if int(len(goldenList.intersection(geneList))/len(goldenList)*100) >= percentage:
        return True
    return False 

def getEmptyDict(rawGoldenList):
    rawGoldenDict = {gene:0 for gene in rawGoldenList}
    return rawGoldenDict

def getGoldenDictValues(rawGoldenList, rawGoldenDict, percentage: int):
    validatedSpecieCount= 0
    failedCount= 0
    for filename in os.listdir('FiltheredGeneLists'):
        listPath= os.path.join('FiltheredGeneLists', filename)
        with open(listPath, 'r') as genesList:
            genesList= [nameList.rstrip() for nameList in genesList]
        setGeneList= set(genesList)
        if not(speciesValidation(genesList, rawGoldenList, percentage)):
            os.remove(listPath)
            failedCount+= 1
            continue
        validatedSpecieCount+= 1
        setGeneList= rawGoldenList.intersection(setGeneList)
        for gene in setGeneList:
            rawGoldenDict[gene] +=1
    return rawGoldenDict, validatedSpecieCount, failedCount

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
    return None

def geneValidation(geneValue, validatedSpecieCount, percentage: int): #8
    if percentage > 100 or percentage < 0:
        sys.exit("error: similarity(%) > 100")
    if geneValue/validatedSpecieCount*100 >= percentage:
        return True
    else: return False

def getFiltheredGoldenDict(goldDict,validatedSpecieCount,geneValidationThreshhold): #8
    for gene in list(goldDict.keys()):
        value = goldDict[gene]
        if(not(geneValidation(value,validatedSpecieCount,geneValidationThreshhold))):
            del goldDict[gene]
    return goldDict  

def verifyDir(directory,unfiltheredDir, pathToGoldList):
    filenames = [filename for filename in os.listdir(directory)]
    pathToGoldGene = pathToGoldList.replace("filteredProximity_GeneLists/",'')
    if len(filenames) == 1 and pathToGoldGene in filenames:
        shutil.rmtree(directory)
        sys.exit(f"Error: Only input term {pathToGoldGene} have this level of proximity(%) and similarity(%)")
    
    #Comentar linha a baixo para não apagar dir com genes filtrados por proximidade    
    #shutil.rmtree(directory)
    #shutil.rmtree(unfiltheredDir)
    return None


# User Args
#term,taxonRank, specieValidationThreshhold,geneValidationThreshhold= getUserArguments()
term,taxonRank, specieValidationThreshhold,geneValidationThreshhold='Canis lupus', 'order', 15, 5
# Copy GeneLists Dir.
directory, unfiltheredDir= getLoadindDir('GeneLists')
# sets golden List
rawGoldenList, goldListPath= setGoldenList(term)
# sets dictionary with {gene: count}
rawGoldenDict= getEmptyDict(rawGoldenList)
# sets empty count values for dict. This will also filther all species that fail in validation.
rawGoldenDict, validatedSpecieCount, failedCount= getGoldenDictValues(rawGoldenList, rawGoldenDict, specieValidationThreshhold)
# create file FilthereGeneLists with the name of all species that passed the 1st validation.
getFiltheredSpeciesNameList('GeneLists')
# sets dictionary of all genes that passed the validation
goldenDict= getFiltheredGoldenDict(rawGoldenDict,validatedSpecieCount,geneValidationThreshhold)
e.printOutputToFile('goldenGenesDictionary', 'w', goldenDict)
# 
verifyDir(directory, unfiltheredDir, goldListPath)

