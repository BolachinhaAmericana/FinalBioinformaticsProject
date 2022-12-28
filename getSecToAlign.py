#!/usr/bin/env python3

#Made by Marine and Silva
from essentials import eSearch, eFetch
import shutil
import os 

def createSquences_Fasta_Dir(directory):
    while True:
        try:
            os.mkdir(directory)
            pathToFastaDir = directory
            return pathToFastaDir
        except FileExistsError:
            shutil.rmtree('Squences_Fasta') 

def getFastas(pathToFastaDir,db):
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







