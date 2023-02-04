#!/usr/bin/env python3
''' Takes one gene at a time and downloads all species with said gene.'''

import shutil
import os
from essentials import entrez_search, entrez_fetch

def fasta_dir_creator(directory):
    '''Created a file for the fastas to go to.'''
    while True:
        try:
            os.mkdir(directory)
            fasta_dir_path = directory
            return fasta_dir_path
        except FileExistsError:
            shutil.rmtree('Squences_Fasta')

def fasta_downloader(fasta_dir_path):
    '''
    takes one gene at a time from FiltredGeneNames list file.
    downloads a fasta file for each gene with the information about said gene
        for all species in FiltheredScientificNames list file.
    '''
    gene_file  = open("FiltredGeneNames_list.txt", "r", encoding='utf8')
    for gene in gene_file:
        names_file = open("FiltredScientificNames_list.txt", "r", encoding='utf8')
        print(gene)
        fasta_path = f'{gene}.fasta'.replace("\n", "")
        for name in names_file:
            while True:
                try:
                    term = f'{name}[Organism] AND {gene}[gene]'
                    search_result = entrez_search('Nucleotide', term)
                    web_env = search_result["WebEnv"]
                    query_key = search_result["QueryKey"]

                    print(f"Downloading record {name} of {gene}")

                    data = entrez_fetch('Nucleotide', web_env,query_key)

                    fasta_file = open(fasta_path, "a", encoding='utf8')
                    fasta_file.writelines(f"{data}\n")
                    fasta_file.close()
                    break
                except Exception:
                    print('Error while downloading. Trying again...')
        shutil.copy2(fasta_path,fasta_dir_path)
        os.remove(fasta_path)
    gene_file.close()
    names_file.close()

if __name__ == '__main__':
    WORKING_PATH = fasta_dir_creator(directory = 'Squences_Fasta')
    fasta_downloader(WORKING_PATH)
