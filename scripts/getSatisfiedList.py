'''filthers the list of genes by species and genes.'''
#!/usr/bin/env python3
# NOW
import shutil
import sys
import os
from essentials import get_user_arguments


def directory_loader(path_gene_lists_dir):
    '''
    What for:
        Create a new dir, and copy all the files from the 'path_gene_lists_dir' into this new dir.
        If the 'filteredProximity_GeneLists' already exists,
        it will delete the directory and try creating it again.
    Arguments:
        path_gene_lists_dir: The path to the directory containing the gene lists to be copied.

    Vars:
        directory: string representing the name of the newly created directory.
        filename: representing the name of the file being copied from the 'path_gene_lists_dir'.

    Returns:
        directory: The name of the newly created directory.

    '''
    while True:
        try:
            os.mkdir('filteredProximity_GeneLists')
            directory = 'filteredProximity_GeneLists'
            for filename in os.listdir(path_gene_lists_dir):
                shutil.copy2(os.path.join(path_gene_lists_dir,filename), directory)
            return directory
        except FileExistsError:
            shutil.rmtree('filteredProximity_GeneLists')


def target_list_finder(search_term,directory):
    '''
    What for:
        Get the gene list for the specified search_term.

        search_term: The search_term for which to get the gene list.
        directory: The directory containing the gene lists.
    Vars:
        path_list: This string represents the path to the file being processed in the for loop.
        path_target_list: This string represents the path to the gene list file.
        target_file_list_genes_list:
            This is a list of strings representing the genes names from the gene list file.
        target_list: This is a set of the target_file_list_genes_list,
            representing the gene names from the gene list file.
    Returns:
        A tuple containing the set of genes in the gene list (target_list),
            the directory containing the gene lists (directory),
            and the path to the gene list file (path_target_list).
    '''
    try:
        for filename in os.listdir(directory):
            path_list = os.path.join(directory, filename)
            path_target_list = f"{directory}/{search_term}_GeneList"
            if path_list == path_target_list:
                with open(path_target_list, 'r', encoding='utf8') as target_file_list:
                    target_file_list_genes_list = [nameList.rstrip() for nameList in target_file_list]
                target_list = set(target_file_list_genes_list)
        return target_list, directory ,path_target_list
    except UnboundLocalError:
        print("Bad Specie name")

def proximity_tester(directory, species_genes_list, target_list, proximity_value):
    '''
    What for:
        Calculates the proximity_value percentage of the Gold gene
        related to the genes of other organisms of the list
    Arguments:
        directory: wanted directory
        species_genes_list: List of genes
        target_list: Animal of input
        proximity_value: proximity_value percentage
    Vars:
        none
    Returns:
        logical value if the proximity_value percentage calculated is
        equal or superior of the proximity_value percentage inputed
    '''
    proximity_value = int(proximity_value)
    if proximity_value > 100 or proximity_value < 0:
        shutil.rmtree(directory)
        sys.exit("error: proximity_value(%) > 100")
    if int(len(target_list.intersection(species_genes_list))/len(target_list)*100) >= proximity_value:
        return True
    return False

def set_target_genes_empty_dict(target_list):
    '''
    What for:
        Create an empty dictionary
            with the genes in 'target_list' as keys and values initialized to 0.
    Arguments:
        target_list: The set of genes to use as keys in the dictionary.
    Vars:
        target_genes_dict:
            A dictionary with the genes in 'target_list' as keys and values initialized to 0.
    Returns:
        target_genes_dict
    '''
    target_genes_dict = {gene:0 for gene in target_list}
    return target_genes_dict

def set_target_dict_values(target_list, directory, target_genes_dict, proximity_value):
    '''
    What for:
        Update the values in the 'target_genes_dict' dictionary
            for the genes that intersect with the gene lists in the 'directory'.
    Arguments:
        target_list: Set of genes in the gold gene list.
        directory: Directory containing the gene lists.
        target_genes_dict: Dictionary with the genes in 'target_list'
            as keys and values to be updated.
        proximity_value: proximity_value value.
    Vars:
        intersected_count:
            Integer that keeps track of the number of gene lists that intersect with the gold list.
        path_list: String that represents the file path to a file in the directory.
        gene_list: List that represents the genes in a gene list file.
        species_genes_list: Represents the gene list read from a file in the directory.
        gene: Represents a gene in the intersection of the gold list
            and the gene list read from a file in the directory.
    Returns:
        A tuple containing the updated 'target_genes_dict' dictionary
            and the count of intersecting gene lists.
    '''
    intersected_count = 0
    for filename in os.listdir(directory):
        path_list = os.path.join(directory, filename)
        with open(path_list, 'r', encoding = 'utf8') as gene_list:
            gene_list = [nameList.rstrip() for nameList in gene_list]
        species_genes_list = set(gene_list)
        if not proximity_tester(directory,species_genes_list, target_list,proximity_value):
            os.remove(path_list)
            continue
        intersected_count += 1
        species_genes_list = target_list.intersection(species_genes_list)
        for gene in species_genes_list:
            target_genes_dict[gene] += 1
    return  target_genes_dict, intersected_count

def get_approved_scientific_names_list(directory):
    '''
    What for:
        Read the names of all the files in a given directory and rewrite them
            in 'FiltredScientificNames_list.txt'.
        The names of the files are assumed to be scientific names.

    Arguments:
        directory: The path to the directory where the scientific names are stored.

    Vars:
        path_list: Represents the file path to a file in the directory.
        scientific_name: Represents the scientific name of an organism, extracted from the file.
        filtredScientificNames_list:
            Represents the file where the scientific names are being written.
    Returns:
        none
    '''
    while True:
        try:
            #os.remove("ScientificNames_list.txt")
            for filename in os.listdir(directory):
                path_list = os.path.join(directory, filename)
                path_list = path_list.replace("filteredProximity_GeneLists/",'')
                scientific_name = path_list.replace("_GeneList",'')
                filthered_scientific_names_list = open("FiltredScientificNames_list.txt", "a", encoding='utf8')
                filthered_scientific_names_list.writelines(f"{scientific_name}\n")
                filthered_scientific_names_list.close()
            break
        except FileNotFoundError:
            open('ScientificNames_list.txt', "x", encoding='utf8')
        except FileExistsError:
            os.remove("FiltredScientificNames_list.txt")

def similarity_tester(directory,value,intersected_count,similarity_value):
    '''
    What for:
        Determine whether the value of a given metric is above a specified threshold.
    Arguments:
        directory: Path to a directory.
        value: Value of the metric to be compared to the threshold.
        intersected_count: Number of elements in the set being compared.
        similarity_value: Threshold value to be compared to the metric value.
    Vars:
        none
    Returns:
        True if the value of the metric is above the threshold, False otherwise.
    '''
    similarity_value = int(similarity_value)
    if similarity_value > 100 or similarity_value < 0:
        shutil.rmtree(directory)
        sys.exit("error: similarity_value(%) > 100")
    if value/intersected_count*100 >= similarity_value:
        return True
    return False

def get_approved_genes_list(directory,target_genes_dict,intersected_count,similarity_value):
    '''
    What for:
        Write the names of genes that pass a given threshold to a file.
    Arguments:
        directory: Path to a directory.
        value: Value of the metric to be compared to the threshold.
        target_genes_dict:
            Dictionary that stores gene names as keys and the number of occurrences as values.
        intersected_count: Number of elements in the set being compared.
        similarity_value: Threshold value to be compared to the metric value.
    Vars:
        value: Rrepresents the number of occurrences of a gene in the target_genes_dict dictionary.
        FiltredGeneNames_list:
            Represents the file where the gene names that pass the threshold are being written.
    Returns:
        target_genes_dict: Updated dictionary.
    '''
    while True:
        try:
            os.remove("FiltredGeneNames_list.txt")
            for gene in list(target_genes_dict.keys()):
                value = target_genes_dict[gene]
                if not similarity_tester(directory,value,intersected_count,similarity_value):
                    del target_genes_dict[gene]
                    continue
                get_approved_gene_list = open("FiltredGeneNames_list.txt", "a", encoding='utf8')
                get_approved_gene_list.writelines(f"{gene}\n")
                get_approved_gene_list.close()
            break
        except FileNotFoundError:
            open('FiltredGeneNames_list.txt', "x", encoding='utf8')
    return target_genes_dict

def dir_checker(directory,path_target_list):
    '''
    What for:
        Check if a directory contains only the gold list file. If it does
            delete the directory and exit the program. Otherwise, delete the directory.
    Arguments:
        directory: Path to the directory to be checked.
        path_target_list: File path to the gold list.
    Vars:
        path_target_genes:
            Represents the name of the gold list file, extracted from path_target_list.
    Returns:
        none
    '''
    filenames = [filename for filename in os.listdir(directory)]
    path_target_genes = path_target_list.replace("filteredProximity_GeneLists/",'')
    if len(filenames) == 1 and path_target_genes in filenames:
        shutil.rmtree(directory)
        sys.exit(f"Error: Only input search_term {path_target_genes} have this level of proximity_value(%) and similarity_value(%)")
    #Comment the line below to not delete the directory with filtered genes by proximity_value.
    shutil.rmtree(directory)

if __name__ == "__main__":
    DIRECTORY = directory_loader(path_gene_lists_dir='GeneLists')
    term,proximity,similar_value = get_user_arguments(3)
    golden_list, folder, path_gold_list = target_list_finder(term, DIRECTORY)
    golden_dict = set_target_genes_empty_dict(golden_list)
    golden_dict, intersect_count = set_target_dict_values(golden_list, folder, golden_dict, proximity)
    get_approved_scientific_names_list(folder)
    filtredDict = get_approved_genes_list(folder, golden_dict, intersect_count, similar_value)
    dir_checker(folder,path_gold_list)
