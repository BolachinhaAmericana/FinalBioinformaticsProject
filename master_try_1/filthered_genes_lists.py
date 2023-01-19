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
def get_user_arguments(): #
    '''
    getUserArguments(term,proximity,similarity):
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
    getsetGoldList(term,directory)
    return setGoldList, directory ,pathToGoldList
    '''
    try:
        for file_name in os.listdir(working_directory):
            path_list = os.path.join(working_directory, file_name)
            path_target_list = f"{working_directory}/{search_term}_GeneList"

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
    intersect(directory,setGeneList,setGoldList,proximity):
    return T/F
    '''
    proximity_value= int(proximity_value)
    if int(len(target_list.intersection(species_genes_list))/len(target_list)*100)>=proximity_value:
        return True
    return False

def set_target_genes_empty_dict(target_list):
    '''
    getGoldEmptyDict(setGoldList)
    return goldDict
    '''
    target_genes_dict = {gene:0 for gene in target_list}
    return target_genes_dict

def set_target_dict_values(target_list, working_directory, target_genes_dict, proximity_value: int):
    '''
    getGoldDictValue(setGoldList, directory, goldDict, proximity):
    return goldDict, intresectCount
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
    getFilteredScientificName_list(directory)
    returns FiltheredScientificnames_list.txt
    '''
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
    haveSimilarity(directory,value,intresectCount,similarity)
    return true/false
    '''
    similarity_value= int(similarity_value)
    if dictionary_value/intersected_count*100 >= similarity_value:
        return True
    return False

def get_approved_genes_list(working_directory,target_genes_dict,intersected_count,similarity_value):
    '''
    FiltredGeneNames_list(directory,goldDict,intresectCount,similarity)
    return goldDict
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
    verifyDir(directory,pathToGoldList)

    '''
    filenames = [filename for filename in os.listdir(working_directory)]
    path_target_genes = path_target_list.replace("filteredProximity_GeneLists/",'')
    if len(filenames) == 1 and path_target_genes in filenames:
        shutil.rmtree(working_directory)
        sys.exit(f"Error: Only input search_term {path_target_genes} have this level of proximity(%) and similarity(%)")


if __name__ == "__main__":

    directory_loader('./GeneLists')

    search,proximity,similarity = get_user_arguments()

    setGoldList, pathToGoldList = target_list_finder(search,working_dir)
    goldDict = set_target_genes_empty_dict(setGoldList)
    goldDict,intresectCount = set_target_dict_values(setGoldList, working_dir, goldDict,proximity)
    get_approved_scientific_names_list(working_dir)
    filtredDict = get_approved_genes_list(working_dir,goldDict,intresectCount,similarity)
    #print(filtredDict)
    dir_checker(working_dir,pathToGoldList)