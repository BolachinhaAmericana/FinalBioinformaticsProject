#!/usr/bin/env python3

#Made by Marine and Silva
from essentials import eSearch, eFetch
import shutil
import os 

def createSquences_Fasta_Dir(directory):
    '''
    What for:
        Creates a directory at the specified path.
        If the directory already exists, it will be removed and a new one will be created.

    Arguments:
        directory: The path to the directory to be created

    Vars:
        pathToFastaDir: Stores the path to the created directory

    Returns:
        pathToFastaDir: The path to the created directory

    '''
    while True:
        try:
            os.mkdir(directory)
            pathToFastaDir = directory
            return pathToFastaDir
        except FileExistsError:
            shutil.rmtree('Squences_Fasta') 

def getFastas(pathToFastaDir,db):
    '''
    What for:
        Takes a path to a directory and a database name as inputs, and retrieves fasta sequences for a list of genes and scientific names.
        This function also reads two files, 'FiltredGeneNames_list.txt' and 'FiltredScientificNames_list.txt' and uses the eSearch and eFetch to fetch the fasta sequences, and save them into a file in the wanted path.

    Arguments:
        pathToFastaDir: The path to the directory where the fasta files will be saved.
        db: The name of the database to search for the fasta sequences.

    Vars:
        gene_file: file that contains the list of gene names.
        scientificNames_file: file that contains the list of scientific names.
        fastaPath: Contains the path to the fasta file.
        term: Contains the search term to fetch the fasta sequence, it is constructed by concatenating scientificName and gene
        eSearchResults: Contains the result of the eSearch function, it is used to fetch the fasta sequence
        webEnv: Contains the web environment needed to fetch the fasta sequence
        queryKey: Contains the query key needed to fetch the fasta sequence
        data: Contains the fasta sequence 
        fastaFile: File used to write the fasta sequence

    Returns:
        none

    '''
    gene_file  = open(f"FiltredGeneNames_list.txt", "r")
    for gene in gene_file:
        scientificNames_file = open("FiltredScientificNames_list.txt", "r")
        fastaPath = f'{gene}.fasta'
        print(gene)
        for scientificName in scientificNames_file:
            while True:
                try:
                    term = f'{scientificName}[Organism] AND {gene}[gene]'
                    eSearchResults = eSearch(db,term)
                    webEnv = eSearchResults["WebEnv"]
                    queryKey = eSearchResults["QueryKey"]
                    print(f"Downloading record {scientificName} of {gene}")
                    data = eFetch(db,webEnv,queryKey)
                    fastaFile = open(fastaPath, "a")
                    fastaFile.writelines("%s\n" %data)
                    fastaFile.close() 
                    break
                except:
                    print('Error while downloading. Trying again...')                       
        shutil.copy2(fastaPath,pathToFastaDir)
        os.remove(fastaPath)        
    gene_file.close()
    scientificNames_file.close()
            


if __name__ == '__main__':
    pathToFastaDir = createSquences_Fasta_Dir(directory = 'Squences_Fasta')
    getFastas(pathToFastaDir,db='Nucleotide')








