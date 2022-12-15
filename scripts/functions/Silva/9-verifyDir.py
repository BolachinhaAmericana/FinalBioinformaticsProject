
import sys
import shutil
import os




def verifyDir(directory,pathToGoldList):
    filenames = [filename for filename in os.listdir(directory)]
    pathToGoldGene = pathToGoldList.replace("filteredProximity_GeneLists/",'')
    if len(filenames) == 1 and pathToGoldGene in filenames:
        shutil.rmtree(directory)
        sys.exit(f"Error: Only input term {pathToGoldGene} have this level of proximity(%) and similarity(%)")
    #Comentar linha a baixo para não apagar dir com genes filtrados por proximidade    
    #shutil.rmtree(directory)