#!/usr/bin/env python3
'''
Will read the list of scientific names and download all gene names from those species.
'''

import time
import concurrent.futures
import os
import shutil
import sys
from Bio import Entrez
from essentials import entrez_search, file_reader_to_list, print_to_file




def species_gene_list_downloader(species: str, batch_size: int):
    '''
    What for:
        This function is used to download a list of genes from the NCBI Gene database for a given species. 
        Uses the Entrez API to search the database and retrieve the genes. 
        To avoid overloading the API, the genes are retrieved in batches using the batch_size parameter.
        Uses multiprocessing while fetching the information.

    Arguments:
        species: Species for which the genes are to be retrieved.
        batch_size: Number of genes to be retrieved at a time.

    Vars:
        search_result: Search results obtained from the Entrez API
        count: Total number of genes found for the species
        web_environment: Web environment parameter required for the Entrez API
        query_key: Query key parameter required for the Entrez API
        gene_list: Contains the names of the genes retrieved from the API

    Returns:
        None

    '''
    search_result = entrez_search('Gene', f'{species} AND alive[prop]')
    count = int(search_result["Count"])
    web_environment = search_result["WebEnv"]
    query_key = search_result["QueryKey"]
    gene_list = set()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(fetch_genes, web_environment, query_key, start, batch_size, species, count) for start in range(0, count, batch_size)]
        for future in concurrent.futures.as_completed(futures):
            gene_list |= future.result()
    try:
        with open(f'./GeneLists/{species}_GeneList.txt', 'w', encoding='utf8') as f:
            for gene in gene_list:
                f.write(gene + '\n')
        return None
    except FileNotFoundError:
        print(f'{species} has 0 genes.')
        print_to_file('./GeneLists/noGene.txt', 'a', species)



def fetch_genes(web_environment, query_key, start, batch_size, species, count):
    '''
    What for:
        This function is used to fetch a batch of gene records from NCBI Gene database for a specific species.
        It takes a batch of gene records starting from the 'start' index, up to 'batch_size' number of records.

    Arguments:
        web_environment: Web environment of the Entrez search results.
        query_key: Query key of the Entrez search results.
        start: Starting index of the batch of records to fetch.
        batch_size: Number of records to fetch in the batch.
        species: Species name to search for in the NCBI Gene database.
        count: Total number of records in the search results.

    Vars:
        genes: A set to store the fetched genes.

    Returns:
       genes: A set of gene names fetched from the NCBI Gene database.
        
    '''
    genes = set()
    while True:
        try:
            print(f"Downloading record {start+1} to {start+batch_size} of {count} - {species}")
            fetch_handle = Entrez.efetch(db="Gene",
                                         rettype="gb",
                                         retmode="text",
                                         retstart=start,
                                         retmax=batch_size,
                                         webenv=web_environment,
                                         query_key=query_key,
                                         idtype="acc")
            print(f"Completed: {start+1} to {start+batch_size} of {count} - {species}")
            time.sleep(1)

            data = fetch_handle.read()
            fetch_handle.close()
            # working data after gathering
            for line in data.splitlines():
                if ']' not in line and ':' not in line:
                    filthered = line.strip().split(' ')[-1]
                    genes.add(filthered)
            
            break
        except Exception:
            print('Error while downloading. Trying again...')
    return genes

SPECIES_NAMES_LIST = file_reader_to_list('./ScientificNames_list.txt')

def downloading():
    '''
    What for:
        This function is used to run the species_gene_list_downloader function for a list of species.

    Arguments:
        None

    Vars:
        SPECIES_NAMES_LIST: List of species names used as input for the species_gene_list_downloader function.
    
    Returns:
        None

    '''
    for specie in SPECIES_NAMES_LIST:
        species_gene_list_downloader(specie, 7500)



def empty_downloads_handler():
    '''
    What for:
        This function is used to handle species that do not have any genes and add their names to the 'noGenes.txt' file.

    Arguments:
        None

    Vars:
        SPECIES_NAMES_LIST: List of species names
        specie: Variable that holds the current species name

    Returns:
        None
        
    '''
    for specie in SPECIES_NAMES_LIST:
        try:
            if os.path.getsize(f'./GeneLists/{specie}_GeneList.txt') == 0:
                os.remove(f'./GeneLists/{specie}_GeneList.txt')
                print_to_file('./GeneLists/noGene.txt', 'a', specie)
                print(f'{specie} has 0 genes.')
        except Exception:
            print('idk why this wouldnt work but yeah...')
            continue

if __name__ == "__main__":
    try:
        if os.path.exists('./GeneLists'):
            shutil.rmtree('./GeneLists')
            os.rmdir('./GeneLists')
    except FileNotFoundError:
        pass
    except OSError:
        print('Please delete GeneLists dir manually')
        sys.exit()
    os.mkdir('./GeneLists')
    downloading()
    empty_downloads_handler()
