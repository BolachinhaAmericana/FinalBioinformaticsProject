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
from essentials import entrez_search, print_to_file, file_reader_to_list

def env_setup(output_dir):
    '''checks for existing and missing files'''
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        shutil.rmtree(output_dir)
        os.rmdir(output_dir)
        os.mkdir(output_dir)
    except FileNotFoundError:
        print('File Not Found. Continuing...')
    except OSError:
        print(f'An Error occured while deleting {output_dir} directory, please remove it manually.')
        sys.exit()
    return None

def species_gene_list_downloader(species: str, batch_size: int):
    '''
    gets all information needed from the Entrez API
    uses multiprocessing while fetching the information.
    '''
    search_result = entrez_search('Gene', f'{species}[Organism] AND alive[prop]')
    count = int(search_result["Count"])
    web_environment = search_result["WebEnv"]
    query_key = search_result["QueryKey"]
    gene_list = set()


    def fetch_genes(web_environment, query_key, start, batch_size, species, count):
        """Function to fetch genes for a specific batch"""
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



    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(fetch_genes, web_environment,
                                                query_key,
                                                start,
                                                batch_size,
                                                species,
                                                count)
                                                for start in range(0, count, batch_size)]

        for future in concurrent.futures.as_completed(futures):
            gene_list |= future.result()
        if '' in gene_list:
            gene_list.remove('')
    try:
        if len(gene_list) == 0:
            print_to_file('./GeneLists/noGene', 'a', species)
            print(f'{species} has 0 genes.')
        else:
            with open(f'./GeneLists/{species}_GeneList', 'w', encoding='utf8') as file:
                for gene in gene_list:
                    file.write(gene + '\n')
        return None
    except FileNotFoundError:
        print('An error occured')

if __name__ == "__main__":
    OUTPUT_DIR = './GeneLists'
    try:
        SPECIES_NAMES_LIST = file_reader_to_list('./ScientificNames_list.txt')
    except FileNotFoundError:
        print('./ScientificNames_list file Not found. Cjeck working dir. Exiting.')
        sys.exit()

    env_setup(OUTPUT_DIR)
    for specie in SPECIES_NAMES_LIST:
        species_gene_list_downloader(specie, 7500)
