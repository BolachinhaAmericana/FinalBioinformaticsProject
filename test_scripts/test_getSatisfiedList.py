import os
import shutil
from getSatisfiedList import getLoadingDir, getsetGoldList, intersect, getGoldEmptyDict

def test_getLoadingDir():
    # Create a test directory with some test files
    test_dir = 'test_dir'
    os.mkdir(test_dir)
    test_files = ['test1.txt', 'test2.txt', 'test3.txt']
    for file in test_files:
        open(os.path.join(test_dir, file), 'w').close()
    
    # Use the function to copy the test files to a new directory
    new_dir = getLoadingDir(test_dir)
    
    # Check that the new directory was created
    assert os.path.exists(new_dir) == True
    
    # Check that the test files were copied to the new directory
    for file in test_files:
        assert os.path.exists(os.path.join(new_dir, file)) == True
    
    # Clean up by removing the test directory and new directory
    shutil.rmtree(test_dir)
    shutil.rmtree(new_dir)

def test_getsetGoldList():
    # Criação de conjunto de ouro de teste
    gold_set = {'gene1', 'gene2', 'gene3'}
    
    # Execução da função
    test_set = getsetGoldList('test_term')
    
    # Verificação se o conjunto retornado é igual ao conjunto de ouro de teste
    assert test_set == gold_set

def test_intersect():
    # Criação de conjuntos de teste
    set1 = {'gene1', 'gene2', 'gene3'}
    set2 = {'gene2', 'gene3', 'gene4'}
    
    # Execução da função
    percentage = intersect(set1, set2)
    
    # Verificação se a porcentagem calculada é correta
    assert percentage == (2/3)*100

def test_getGoldEmptyDict():
    # Criação de conjunto de ouro de teste
    gold_set = {'gene1', 'gene2', 'gene3'}
    
    # Execução da função
    test_dict = getGoldEmptyDict()
    
    # Verificação se o dicionário criado tem as chaves corretas e os valores inicializados como 0
    assert set(test_dict.keys()) == gold_set
    assert all(value == 0 for value in test_dict.values())