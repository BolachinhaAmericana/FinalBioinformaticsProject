#!/usr/bin/env python3
'''
This script will receive 2 inputs that will be one folder
with all the fasta files and the second is the list of species names.
'''

import subprocess
import shutil
import glob
import os
from Bio import SeqIO

def get_named_fastas(directory):
    '''
    What for:
        This function renames the fasta files in the directory using species names
        from the text file "FiltredScientificNames_list.txt" and saves them in the inputed dir.
        If the inputed directory already exists, it will be removed,
        and a new directory will be created.

    Arguments:
        directory: Directory where the renamed fasta files will be saved.

    Vars:
        new_names: List of scientific names read from the file "FiltredScientificNames_list.txt".
        fasta_files: List of strings representing
        the file paths of all fasta files in the directory "Squences_Fasta".
        records: List of SeqIO records for a single fasta file.
        new_file_name: String representing the new name of the fasta file,
        which is based on the original file name with the "_updated".

    Returns:
       none

    '''
    with open("FiltredScientificNames_list.txt", encoding= 'utf8') as doc:
        new_names = doc.read().splitlines()
    if os.path.exists(directory):
        shutil.rmtree(directory)
    if not os.path.exists(directory):
        os.makedirs(directory)
    fasta_files = glob.glob("Squences_Fasta/*.fasta")
    for fasta in fasta_files:
        records = list(SeqIO.parse(fasta, "fasta"))
        for record, new_name in zip(records, new_names):
            record.id = new_name.replace(" ", "_")
            record.id = new_name.replace("'", "")
            record.description = ""
        new_file_name = os.path.basename(fasta).split(".")[0] + "_updated.fasta"
        with open(f"{directory}/{new_file_name}", "w", encoding= 'utf8') as doc:
            SeqIO.write(records, doc, "fasta")

def fastas_aligner(input_folder,out_folder):
    '''
    What for:
        This function aligns the fasta files in a
        specified input folder using the MAFFT alignment tool.

    Arguments/Vars:
        input_folder: Directory containing the fasta files to be aligned.
        out_folder: Directory where the aligned fasta files will be saved.

    Returns:
       none

    '''
    folder = input_folder
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)
    fasta_files = [doc for doc in os.listdir(folder) if doc.endswith('.fasta')]
    for fasta in fasta_files:
        doc = open(f"{out_folder}/_aligned{fasta}","w", encoding='utf8')
        subprocess.run(['mafft', folder + '/' + fasta], stdout = doc)

def fastas_concatenator(input_folder):
    '''
    What for:
        This function concatenates all the fasta files in the directory and stores the
        concatenated fasta in a new file called "concat.fasta" in the current directory.
        The function also deletes the input_folder after concatenation is done.

    Arguments:
        input_folder: Path to a directory containing multiple fasta files.
    Vars:
        folder: Holds the input folder.
        concatenated_fasta: Dictionary that will hold the concatenated fasta sequences.
        seq: Variable that holds the current sequence read from the fasta file.
        name: Variable that holds the current fasta name read from the fasta file
    Returns:
        None
    '''
    folder = input_folder
    concatenated_fasta = {}
    seq = ""
    name = ""
    for fasta in os.listdir(folder):
        with open (f"{folder}/{fasta}","r", encoding='utf8') as read:
            for line in read:
                if line.startswith(">"):
                    name = line.strip()
                    if name not in concatenated_fasta:
                        concatenated_fasta[name] = []
                else:
                    seq = line.strip()
                    concatenated_fasta[name].append(seq)
    subprocess.run(f'rm -r {folder}', shell = True)
    with open ("concat.fasta","w", encoding='utf8') as doc:
        for key, value in concatenated_fasta.items():
            doc.write(f'{key}\n')
            for i in value:
                doc.write(i+"\n")

if __name__ == "__main__":
    get_named_fastas('named_Fastas')
    fastas_aligner('named_Fastas','alignedFastas')
    fastas_concatenator('alignedFastas')
