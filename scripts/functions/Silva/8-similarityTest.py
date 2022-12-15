import shutil
import sys

#def haveSimilarity(directory,value,intresectCount,similarity):
#    similarity = int(similarity)
#    if similarity > 100:
#        shutil.rmtree(directory)
#        sys.exit("error: similarity(%) > 100")  
#    if value/intresectCount*100 >= similarity:
#        return True
#    return False

def geneValidation(geneValue, validatedSpecieCount, percentage: int):
    if percentage > 100 or percentage < 0:
        sys.exit("error: similarity(%) > 100")
    if geneValue/validatedSpecieCount*100 >= percentage:
        return True
    else: return False


def getsimilarityGoldDict(goldDict,intresectCount,similarity):
    for gene in list(goldDict.keys()):
        value = goldDict[gene]
        if(not(geneValidation(value,intresectCount,similarity))):
            del goldDict[gene]
    return goldDict  



validatedSpecieCount= 15
rawGoldenDict= {'ND2': 15, 'ND6': 15, 'ND4L': 15, 'ATP6': 15, 'ND1': 15, 'COX3': 15, 'CYTB': 15, 'ATP8': 15, 'ND5': 15, 'COX2': 15, 'ND3': 15, 'ND4': 15, 'COX1': 15}

goldenDict = getsimilarityGoldDict(rawGoldenDict,validatedSpecieCount,50)
#print(goldenDict)
