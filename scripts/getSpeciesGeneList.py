#!/usr/bin/env python3
'''
Will read the list of scientific names and download all gene names from those species.
'''
import os
from Bio import Entrez
from essentials import entrez_search, print_to_file, file_reader_to_list

def species_gene_list_downloader(species: str, batch_size: int): #####
    '''
    Downloads one species gene list
    Arguments:
        specie: scientific name of a specie to get gene list
        batch_size: max= 10000. The number of genes that it will fetch at a time.
    Vars-
        search_result: result of a eSearch, filthering all replaced or discontinued genes.
        count: number of total genes that will be fetched
        web_environment: Web Environment
        query_key: Query Key
        fetch_handle: eFetch Query.
        data: read of the fetch_handle
    Returns: None, but will create a LoadingGeneList.txt file

    '''
    search_result= entrez_search('Gene', f'{species} AND alive[prop]')
    count = int(search_result["Count"])
    web_environment = search_result["WebEnv"]
    query_key = search_result["QueryKey"]

    for start in range(0, count, batch_size):
        while True:
            try:
                print(f"Downloading record {start+1} to {count} - {specie}")
                fetch_handle = Entrez.efetch(db="Gene",
                                             rettype="gb",
                                             retmode="text",
                                             retstart=start,
                                             retmax=batch_size,
                                             webenv=web_environment,
                                             query_key=query_key,
                                             idtype="acc")
                data = fetch_handle.read()
                fetch_handle.close()
                print_to_file('loading.dump', 'w', data)
                os.system("grep -v ':' loading.dump > loading1.dump")
                os.system("grep '\.' loading1.dump > loading2.dump")
                os.system("grep -v '\[' loading2.dump > loading3.dump")
                os.system("sed 's/^.*\ //' loading3.dump >> LoadingGeneList.txt")
                os.system('rm ./*.dump')
                break
            except Exception:
                print('Error while downloading. Trying again...')
    return None #Outputs the file  LoadingGeneList.txt

def species_gene_list_handler(species: str, gene_list_path):
    '''
    Organizes the gene list and reacts incase list is empty.
    Arguments:
        takes the specie name and output path as arguments
    Vars-
        None
    Returns:
        None, but will rename and move a file to GeneLists/
    '''
    try:
        os.rename(gene_list_path, f'./GeneLists/{species}_GeneList')
    except FileNotFoundError:
        print_to_file('./GeneLists/noGeneSpecies.txt', 'a', species)
        print(f'The following species has zero genes: {species}')
    return None
#3- Done- Genes Lists





if __name__ == "__main__":
    speciesList= file_reader_to_list('./ScientificNames_list.txt')


    if os.path.exists('./GeneLists'):
        os.system('rm ./GeneLists/*')
    else: os.system('mkdir ./GeneLists')

    for specie in speciesList:
        species_gene_list_downloader(specie, 7500)
        species_gene_list_handler(specie, './LoadingGeneList.txt')
