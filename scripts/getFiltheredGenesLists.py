#!/usr/bin/env python3
'''
This file will remove all genes that don't make the cut for the analisys
Will create filtheredProximity_geneLists with the list of genes that got approved
Will create FiltheredFeneNames_list with the names of the genes that got approved
Will create FiltheredScientificNames_list with the names of the species that got approved
'''

import os
import shutil
import sys


working_dir = './filteredProximity_GeneLists'

def directory_loader(path_gene_lists_dir): #
    '''
    Use:
        Create a new dir and copy all files from 'geneLists_directory' into this new dir.
        Incase 'filteredProximity_GeneLists' exists, will remove and create dir again.
    Arguments:
        geneLists_directory: The path to the directory containing the gene lists to be copied.
    Vars:
        directory: name of the newly created directory.
        filename: name of the file being copied from the 'geneLists_directory'.
    Returns:
        directory: The name of the newly created directory.
    '''
    try:
        os.mkdir('./filteredProximity_GeneLists')
        working_directory = './filteredProximity_GeneLists'
        for filename in os.listdir(path_gene_lists_dir):
            shutil.copy2(os.path.join(path_gene_lists_dir,filename), working_directory)

    except FileExistsError:
        shutil.rmtree('filteredProximity_GeneLists')
        os.mkdir('./filteredProximity_GeneLists')
        working_directory = './filteredProximity_GeneLists'
        for filename in os.listdir(path_gene_lists_dir):
            shutil.copy2(os.path.join(path_gene_lists_dir,filename), working_directory)


def get_user_arguments():
    '''
    What for:
        Get the user wanted search_term, proximity, and similarity values.

    Arguments/Vars:
        search_term: The user-specified search_term.
        proximity_value: The user-specified proximity.
        similarity_value: The user-specified similarity.

    Returns:
        tuple: A tuple containing the search_term, proximity, and similarity values.
    '''
    try:
        search_term = sys.argv[1]
        proximity_value = int(sys.argv[2])
        similarity_value = int(sys.argv[3])
        if all(i in range(0, 101) for i in (proximity_value, similarity_value)):
            return search_term, proximity_value, similarity_value
        else: raise IndexError
    except IndexError:
        sys.exit('''
        This program Takes 3 arguments:
            search term - Original species scientific name
            Proximity value - percentage(%) value (0-100)
            Similatiry value - percentage(%) value (0-100)
            ### Please note that the percentage values must be int !!!
        ''')

def target_list_finder(search_term,working_directory):
    '''
    What for:
        Get the gene list for the specified search_term.
    Arguments:
        search_term: The search_term for which to get the gene list.
        working_directory: The directory containing the gene lists.
    Returns:
        target_list, path_target_list'''
    try:
        for file_name in os.listdir(working_directory):
            path_list = os.path.join(working_directory, file_name)
            path_target_list = f"{working_directory}/{search_term}_GeneList.txt"

            if path_list == path_target_list:

                with open(path_target_list, 'r', encoding='utf8') as target_file_list:
                    target_file_list_genes_list = [genes_list.rstrip() for genes_list in target_file_list]
                target_list = set(target_file_list_genes_list)

        return target_list, path_target_list
    except UnboundLocalError:
        print("Bad Specie name")
        target_list, path_target_list= None, None
        return target_list, path_target_list

def proximity_tester(species_genes_list,target_list,proximity_value: int):
    '''
    What for:
        Calculates the proximity percentage of the target gene
        related to the genes of other organisms on the list
    Arguments:
        species_genes_list: List of genes to verify
        target_list: input species gene list
        proximity_value: proximity percentage to verify
    Returns:
        boolean value if proximity high enough or not
    '''
    proximity_value= int(proximity_value)
    if int(len(target_list.intersection(species_genes_list))/len(target_list)*100)>=proximity_value:
        return True
    return False

def set_target_genes_empty_dict(target_list):
    '''
    Arguments:
        target_list: The set of genes to use as keys in the dictionary.
    Returns:
        target_genes_dict
        - Dictionary containing genes in target_list as keys and values initialized to 0
    '''
    target_genes_dict = {gene:0 for gene in target_list}
    return target_genes_dict

def set_target_dict_values(target_list, working_directory, target_genes_dict, proximity_value: int):
    '''
    What for:
        Update the values in the target dictory
        dictionary for the genes that intersect with the gene lists in the 'directory'.
    Arguments:
        target_list
        working_directory
        target_genes_dictionary
        proximity_value
    Returns:
        target_genes_dict, intersected count
        target -> Genes that got approved
        intersected count -> number of approved genes.
    '''
    intersected_count = 0
    for filename in os.listdir(working_directory):
        list_path = os.path.join(working_directory, filename)
        with open(list_path, 'r', encoding= 'utf8') as gene_list:
            gene_list = [name_list.rstrip() for name_list in gene_list]
        species_genes_list = set(gene_list)
        if not proximity_tester(species_genes_list,target_list,proximity_value):
            os.remove(list_path)
            continue
        intersected_count += 1
        species_genes_list = target_list.intersection(species_genes_list)
        for gene in species_genes_list:
            target_genes_dict[gene] += 1
    return  target_genes_dict, intersected_count

def get_approved_scientific_names_list(working_directory):
    '''
    What for:
        Read the names of the files in a given dir and rewrite them in 'working_directory'.
        The names of the files are assumed to be scientific names.
    Arguments:
        working_directory: The path to the directory where the scientific names are stored.
    Vars:
        list_path: Represents the file path to a file in the directory.
        scientific_name: Represents the scientific name of an organism, extracted from the file.
        approved_scientific_names_list: file where the scientific names are being written
    Returns: none'''
    while True:
        try:
            for filename in os.listdir(working_directory):
                list_path = os.path.join(working_directory, filename)
                list_path = list_path.replace("filteredProximity_GeneLists/",'')
                scientific_name = list_path.replace("_GeneList",'')
                scientific_name = scientific_name.replace('./', '')
                approved_scientific_names_list = open("FiltredScientificNames_list.txt", "a", encoding='utf8')
                approved_scientific_names_list.writelines(f"{scientific_name}\n")
                approved_scientific_names_list.close()
            break

        except FileNotFoundError:
            open('ScientificNames_list.txt', "x", encoding='utf8')
        except FileExistsError:
            os.remove("FiltredScientificNames_list.txt")

def similarity_tester(dictionary_value,intersected_count,similarity_value: int):
    '''
    What for:
        check if similarity_value gets trough the threshold.
    Arguments:
        dictionary_value: Value of the metric to be compared to the threshold.
        intresect_count: Number of elements in the set being compared.
        similarity_value: Threshold value to be compared to the metric value.
    Returns:
        Boolean Value
    '''
    similarity_value= int(similarity_value)
    if dictionary_value/intersected_count*100 >= similarity_value:
        return True
    return False

def get_approved_genes_list(target_genes_dict,intersected_count,similarity_value):
    '''
    What for:
        Write the names of genes that pass a given threshold to a file.
    Returns:
        target_genes_dict: Updated dictionary.
    '''
    while True:
        try:
            os.remove("FiltredGeneNames_list.txt")
            for gene in list(target_genes_dict.keys()):
                value = target_genes_dict[gene]
                if not similarity_tester(value,intersected_count,similarity_value):
                    del target_genes_dict[gene]
                    continue
                approved_genes_list = open("FiltredGeneNames_list.txt", "a", encoding='utf8')
                approved_genes_list.writelines(f"{gene}\n")
                approved_genes_list.close()
            break
        except FileNotFoundError:
            open('FiltredGeneNames_list.txt', "x", encoding='utf8')
    return target_genes_dict

def dir_checker(working_directory,path_target_list):
    '''
    What for:
        Check if a directory contains only the gold list file.
        If it does, delete the directory and exit the program.
        Otherwise, delete the directory.
    Arguments:
        working_directory: Path to the directory to be checked.
        path_target_list: File path to the gold list.
    Vars:
        path_target_genes: Represents the name of the gold list file, extracted from pathToGoldList.
    Returns:
        none'''
    filenames = [filename for filename in os.listdir(working_directory)]
    path_target_genes = path_target_list.replace("filteredProximity_GeneLists/",'')
    if len(filenames) == 1 and path_target_genes in filenames:
        shutil.rmtree(working_directory)
        sys.exit(f'Error: Only input search_term {path_target_genes} got approved!!!')
if __name__ == "__main__":

    directory_loader('./GeneLists')

    search,proximity,similarity = get_user_arguments()

    setGoldList, pathToGoldList = target_list_finder(search,working_dir)
    goldDict = set_target_genes_empty_dict(setGoldList)
    goldDict,intresectCount = set_target_dict_values(setGoldList, working_dir, goldDict,proximity)
    get_approved_scientific_names_list(working_dir)
    filtredDict = get_approved_genes_list(goldDict,intresectCount,similarity)
    #print(filtredDict)
    dir_checker(working_dir,pathToGoldList)
