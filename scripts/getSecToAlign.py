#!/usr/bin/env python3

#Made by Marine and Silva
from essentials import entrez_search, entrez_fetch
import shutil
import os 

def createSquences_Fasta_Dir(directory):
    '''
    What for:
        This function is used to create a directory to store FASTA files.

    Arguments/Vars:
        directory: Name of the directory to be created.

    Returns:
        pathToFastaDir: Path of the created directory. 
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
        This function downloads FASTA sequences from a given database using gene names and scientific names.
         
    Arguments:
        pathToFastaDir: Directory where the FASTA files will be stored
        db: Database to download the sequences from.
    
    Vars:
        gene_file: File object representing the file containing the gene names.
        scientificNames_file: File object representing the file containing the scientific names.
        term: Represents the search term to be used in the entrez_search function.
        eSearchResults: Dictionary containing the results of the entrez_search function.
        webEnv: Represents the web environment of the search results.
        queryKey: Represents the query key of the search results.
        data: Data returned by the entrez_fetch function.
        fastaPath: Represents the path of the FASTA.
        fastaFile: Represents the FASTA file.

    Returns:

    '''
    gene_file  = open(f"FiltredGeneNames_list.txt", "r")
    for gene in gene_file:
        scientificNames_file = open("FiltredScientificNames_list.txt", "r")
        print(gene)
        fastaPath = f'{gene}.fasta'.replace("\n", "")
        for scientificName in scientificNames_file:
            while True:
                try:
                    term = f'{scientificName}[Organism] AND {gene}[gene]'
                    eSearchResults = entrez_search(db,term)
                    webEnv = eSearchResults["WebEnv"]
                    queryKey = eSearchResults["QueryKey"]
                    print(f"Downloading record {scientificName} of {gene}")
                    data = entrez_fetch(db,webEnv,queryKey)
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








