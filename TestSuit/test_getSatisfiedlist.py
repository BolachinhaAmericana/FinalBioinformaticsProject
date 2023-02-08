import os
import shutil

from essentials import get_user_arguments
from scripts.getSatisfiedList import (directory_loader, target_list_finder, proximity_tester, 
                      set_target_genes_empty_dict)

def test_directory_loader():
    path_gene_lists_dir = "test_dir"
    os.mkdir(path_gene_lists_dir)
    filename = "test_file.txt"
    with open(os.path.join(path_gene_lists_dir, filename), "w") as f:
        f.write("test")
    directory = directory_loader(path_gene_lists_dir)
    assert os.path.exists(directory)
    assert os.path.exists(os.path.join(directory, filename))
    shutil.rmtree(path_gene_lists_dir)
    shutil.rmtree(directory)

def test_target_list_finder():
    search_term = "test_term"
    directory = "test_dir"
    os.mkdir(directory)
    filename = f"{search_term}_GeneList"
    with open(os.path.join(directory, filename), "w") as f:
        f.write("gene1\ngene2")
    target_list, directory, path_target_list = target_list_finder(search_term, directory)
    assert target_list == {"gene1", "gene2"}
    assert path_target_list == os.path.join(directory, filename)
    shutil.rmtree(directory)

def test_proximity_tester():
    directory = "test_dir"
    os.mkdir(directory)
    species_genes_list = {"gene1", "gene2"}
    target_list = {"gene2", "gene3"}
    proximity_value = 50
    result = proximity_tester(directory, species_genes_list, target_list, proximity_value)
    assert result == True
    proximity_value = 75
    result = proximity_tester(directory, species_genes_list, target_list, proximity_value)
    assert result == False
    shutil.rmtree(directory)

def test_set_target_genes_empty_dict():
    target_list = {"gene1", "gene2"}
    target_genes_dict = set_target_genes_empty_dict(target_list)
    assert target_genes_dict == {"gene1": 0, "gene2": 0}


def test_proximity_tester():
    # Test 1: proximity_value = 50, intersection = target_list, expected result = True
    target_list = set(['gene1', 'gene2', 'gene3'])
    species_genes_list = set(['gene1', 'gene2', 'gene3'])
    proximity_value = 50
    directory = 'filteredProximity_GeneLists'
    result = proximity_tester(directory, species_genes_list, target_list, proximity_value)
    assert result == True, f"For proximity_value = {proximity_value}, species_genes_list = {species_genes_list}, target_list = {target_list}, expected result = True but got {result}"

    # Test 2: proximity_value = 75, intersection = target_list, expected result = True
    target_list = set(['gene1', 'gene2', 'gene3'])
    species_genes_list = set(['gene1', 'gene2', 'gene3'])
    proximity_value = 75
    directory = 'filteredProximity_GeneLists'
    result = proximity_tester(directory, species_genes_list, target_list, proximity_value)
    assert result == True, f"For proximity_value = {proximity_value}, species_genes_list = {species_genes_list}, target_list = {target_list}, expected result = True but got {result}"

    # Test 3: proximity_value = 75, intersection = 2, expected result = False
    target_list = set(['gene1', 'gene2', 'gene3'])
    species_genes_list = set(['gene1', 'gene2'])
    proximity_value = 75
    directory = 'filteredProximity_GeneLists'
    result = proximity_tester(directory, species_genes_list, target_list, proximity_value)
    assert result == False, f"For proximity_value = {proximity_value}, species_genes_list = {species_genes_list}, target_list = {target_list}, expected result = False but got {result}"