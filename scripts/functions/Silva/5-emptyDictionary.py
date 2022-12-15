def getEmptyDict(rawGoldenList):
    rawGoldenDict = {gene:0 for gene in rawGoldenList}
    return rawGoldenDict

rawGoldenList= {'ATP8', 'ND1', 'ND3', 'COX3', 'ND5', 'COX1', 'ND6', 'COX2', 'ND2', 'ND4L', 'CYTB', 'ATP6', 'ND4'}
rawGoldenDict= getEmptyDict(rawGoldenList)
print(rawGoldenDict)

