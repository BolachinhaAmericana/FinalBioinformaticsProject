import os 

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

            

term= 'Canis latrans'
#term= 'Ailurus fulgens'
rawGoldenList= setGoldenList(term)
print(rawGoldenList)

