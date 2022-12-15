import sys
import os

def speciesTestOne(geneList, goldenList, percentage: int): #Filthering species.
    if percentage > 100 or percentage < 0:
            sys.exit("error: proximity(%) > 100")
    if int(len(goldenList.intersection(geneList))/len(goldenList)*100) >= percentage:
        return True
    return False 



def getGoldenDictValues(rawGoldenList, rawGoldenDict, percentage: int):
    validatedSpecieCount= 0
    failedCount= 0
    for filename in os.listdir('FiltheredLists'):
        listPath= os.path.join('FiltheredLists', filename)
        with open(listPath, 'r') as genesList:
            genesList= [nameList.rstrip() for nameList in genesList]
        setGeneList= set(genesList)
        if not(speciesTestOne(genesList, rawGoldenList, percentage)):
            os.remove(listPath)
            failedCount+= 1
            continue
        validatedSpecieCount+= 1
        setGeneList= rawGoldenList.intersection(setGeneList)
        for gene in setGeneList:
            rawGoldenDict[gene] +=1
    return rawGoldenDict, validatedSpecieCount, failedCount


rawGoldenList= {'ND4L', 'ND1', 'COX1', 'ATP6', 'COX3', 'COX2', 'ND5', 'ATP8', 'CYTB', 'ND6', 'ND4', 'ND2', 'ND3'}
rawGoldenDict= {'ND2': 0, 'ND6': 0, 'ND4L': 0, 'ATP6': 0, 'ND1': 0, 'COX3': 0, 'CYTB': 0, 'ATP8': 0, 'ND5': 0, 'COX2': 0, 'ND3': 0, 'ND4': 0, 'COX1': 0}
percentage= 20

rawGoldenDict, validatedSpecieCount, failedCount= getGoldenDictValues(rawGoldenList, rawGoldenDict, percentage)
print(rawGoldenDict, validatedSpecieCount, failedCount)

