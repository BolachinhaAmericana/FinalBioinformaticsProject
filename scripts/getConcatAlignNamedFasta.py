#!/usr/bin/env python3
from Bio import SeqIO
import subprocess
import shutil
import glob
import os
# This script will receive 2 inputs that will be one folder with all the fasta files and the second is the list of species names.

# Made by Pinto feat Silva (getNamed_fastas)
def getNamed_fastas(dir):
    '''
    What for:
        This function renames the fasta files in the directory using species names from the text file "FiltredScientificNames_list.txt" and saves them in the inputed directory.
        If the inputed directory already exists, it will be removed, and a new directory will be created.

    Arguments:
        dir: Directory where the renamed fasta files will be saved.
    
    Vars:
        new_names: List of scientific names read from the file "FiltredScientificNames_list.txt".
        fasta_files: List of strings representing the file paths of all fasta files in the directory "Squences_Fasta".
        records: List of SeqIO records for a single fasta file.
        new_file_name: String representing the new name of the fasta file, which is based on the original file name with the "_updated".
        
    Returns:
       none

    '''
    with open("FiltredScientificNames_list.txt") as doc:
        new_names = doc.read().splitlines()
    if os.path.exists(dir):
        shutil.rmtree(dir)
    if not os.path.exists(dir):
        os.makedirs(dir)
    fasta_files = glob.glob("Squences_Fasta/*.fasta")
    for fasta in fasta_files:
        records = list(SeqIO.parse(fasta, "fasta"))
        for record, new_name in zip(records, new_names):
            record.id = new_name.replace(" ", "_") 
            record.description = ""
        new_file_name = os.path.basename(fasta).split(".")[0] + "_updated.fasta"
        with open(f"{dir}/{new_file_name}", "w") as doc:
            SeqIO.write(records, doc, "fasta")



def fastas_Aligner(input_Folder,out_folder):
    '''
    What for:
        This function aligns the fasta files in a specified input folder using the MAFFT alignment tool.

    Arguments/Vars:
        inputFolder: Directory containing the fasta files to be aligned.
        outfolder: Directory where the aligned fasta files will be saved.
        
    Returns:
       none

    '''
    folder = input_Folder
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)
    fasta_Files = [doc for doc in os.listdir(folder) if doc.endswith('.fasta')]
    for fasta in fasta_Files:
        doc = open(f"{out_folder}/_aligned{fasta}","w")
        subprocess.run(['mafft', folder + '/' + fasta], stdout = doc)


def fastas_Concatenator(input_Folder):
    '''
    What for:
        This function concatenates all the fasta files in the directory and stores the concatenated fasta in a new file called "concat.fasta" in the current directory. 
        The function also deletes the inputFolder after concatenation is done.

    Arguments:
        inputFolder: Path to a directory containing multiple fasta files.
    
    Vars:  
        folder: Holds the input folder.
        concatenatedFasta: Dictionary that will hold the concatenated fasta sequences.
        seq: Variable that holds the current sequence read from the fasta file.
        name: Variable that holds the current fasta name read from the fasta file
        
    Returns:
        None

    '''
    folder = input_Folder
    concatenated_Fasta = {}
    seq = ""
    name = ""
    for fasta in os.listdir(folder):
        with open (f"{folder}/{fasta}","r") as read:
            for line in read:
                if line.startswith(">"):
                    name = line.strip()
                    if name not in concatenated_Fasta:
                        concatenated_Fasta[name] = []
                else:
                    seq = line.strip()
                    concatenated_Fasta[name].append(seq)
    subprocess.run(f'rm -r {folder}', shell = True)
    with open ("concat.fasta","w") as doc:
        for key, value in concatenated_Fasta.items():
            doc.write('%s\n' %(key))
            for i in value:
                doc.write(i+"\n")

if __name__ == "__main__":
    getNamed_fastas('named_Fastas')
    fastas_Aligner('named_Fastas','alignedFastas')
    fastas_Concatenator('alignedFastas')
