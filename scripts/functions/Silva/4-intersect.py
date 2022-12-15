import shutil
import sys


def intersect(directory,setGeneList,setGoldList,proximity):
    proximity = int(proximity)
    if proximity > 100:
        shutil.rmtree(directory)
        sys.exit("error: proximity(%) > 100")
    if int(len(setGoldList.intersection(setGeneList))/len(setGoldList)*100) >= proximity:
        return True
    return False 




def speciesTestOne(geneList, goldenList, percentage: int): #Filthering species.
    if percentage > 100 or percentage < 0:
            sys.exit("error: proximity(%) > 100")
    if int(len(goldenList.intersection(geneList))/len(goldenList)*100) >= percentage:
        return True
    return False 

rawGoldenList, goldListPath= {'ATP8', 'ND1', 'ND3', 'COX3', 'ND5', 'COX1', 'ND6', 'COX2', 'ND2', 'ND4L', 'CYTB', 'ATP6', 'ND4'}, 'GeneLists/Canis latrans_GeneList'
geneList= {'COX3', 'trnS(UCN)', 'trnR', 'COX2', 'CYTB', 'BMR74_gp4', 'BMR81_gp2', 'BMR81_gp1', 'trnQ', 'trnW', 'trnL(UUR)', 'ND3', 'trnA', 'rRNA', 'BMR74_gp3', 'trnL(CUN)', 'trnT', 'ATP6', 'trnI', 'ND5', 'BMR74_gp2', 'trnY', 'trnN', 'trnP', 'COX1', 'trnC', 'trnF', 'ND2', 'trnD', 'ND4L', 'trnM', 'trnV', 'BMR74_gp1', 'trnS(AGY)', 'trnK', 'ATP8', 'ND4', 'trnE', 'trnG', 'ND1', 'trnH', 'ND6'}
print(speciesTestOne(rawGoldenList, geneList,100))