import os
import sys
import shutil
from scripts.getSatisfiedList import getLoadindDir, getUserArguments, getsetGoldList, getGoldEmptyDict, getGoldDictValue, getFilteredScientificName_list, FiltredGeneNames_list, verifyDir

#Made by Marine feat Pinto

def test_userArgs():
    sys.argv = ['getUserArguments', '', '', '']
    term, proximity, similarity = getUserArguments(sys.argv[1], sys.argv[2], sys.argv[3])
    assert term == ''
    assert proximity == ''
    assert similarity == ''

def test_getSetGollist():
    pass

def test_getGoldEmpty():
    setGoldList = {'gene1', 'gene2', 'gene3'}
    goldDict = getGoldEmptyDict(setGoldList)
    assert goldDict == {'gene1': 0, 'gene2': 0, 'gene3': 0}

def test_getGolddictValue():
    setGoldList = set(["gene1", "gene2", "gene3"])
    directory = "filteredProximity_GeneLists"
    goldDict = {"gene1": 0, "gene2": 0, "gene3": 0}
    proximity = 10
    result = getGoldDictValue(setGoldList, directory, goldDict, proximity)
    assert result == ({"gene1": 0, "gene2": 0, "gene3": 0}, 0), f"Expected ({{'gene1': 0, 'gene2': 0, 'gene3': 0}}, 0), but got {result}"


def test_getFilteredScientificName_list():
    test_dir = 'test_dir'
    os.mkdir(test_dir)
    open(os.path.join(test_dir, 'file1_GeneList'), 'w').close()
    open(os.path.join(test_dir, 'file2_GeneList'), 'w').close()
    getFilteredScientificName_list(test_dir)
    with open("FiltredScientificNames_list.txt", 'r') as f:
        scientific_names = f.read().strip().split("\n")
    assert type(scientific_names) == list
    assert len(scientific_names) > 0
    shutil.rmtree(test_dir)
    os.remove("FiltredScientificNames_list.txt")   
    
def test_FiltredGeneNames_list():
    directory = "test_directory"
    goldDict = {'gene1': 10, 'gene2': 5, 'gene3': 7}
    intresectCount = 20
    similarity = 0.5
    os.makedirs(directory, exist_ok=True)
    result = FiltredGeneNames_list(directory, goldDict, intresectCount, similarity)
    assert type(result) == dict
    os.rmdir(directory)

def test_verifyDir():
    # create a test directory and a test file
    test_dir = 'test_dir'
    os.mkdir(test_dir)
    open(os.path.join(test_dir, 'test_file.txt'), 'a').close()

    # call the verifyDir function with the test directory and test file
    pathToGoldList = os.path.join(test_dir, 'test_file.txt')
    verifyDir(test_dir, pathToGoldList)

    # check if the test directory has been deleted
    assert not os.path.exists(test_dir)