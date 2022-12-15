import os
import shutil

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

test= getLoadindDir('GeneLists')
directory = getLoadindDir(geneLists_directory='GeneLists')
print(test)