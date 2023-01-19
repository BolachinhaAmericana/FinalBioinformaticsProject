#!/usr/bin/env python3
from Bio import SeqIO
import subprocess
import shutil
import glob
import os

# Made by Pinto feat Silva (getNamed_fastas)


def getNamed_fastas(dir):
    with open("FiltredScientificNames_list.txt") as f:
        new_names = f.read().splitlines()
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
        with open(f"{dir}/{new_file_name}", "w") as f:
            SeqIO.write(records, f, "fasta")

"""
This function has 2 arguments the first one is a folder with the original fastas 
and the second one will be folder where the aligned fastas will be.
If the second folder does not exist, this function will create that directory
and then proceed to align all the fastas
"""
def fastas_Aligner(inputFolder, output_Folder):

    if not os.path.isdir(output_Folder):
        os.makedirs(output_Folder)
    fasta_Files = [f for f in os.listdir(inputFolder)]
    for fasta in fasta_Files:
        f = open(f"{output_Folder}/_aligned{fasta}","w")
        subprocess.run(['mafft', inputFolder + '/' + fasta], stdout = f)


"""
This function will receive the aligned fastas folder and will concatenate them all.
In the end we will get a massive fasta that we will use to build the trees.
"""
def fastas_Concatenator(input_Folder):
    concatenatedFasta = {}
    seq = ""
    seq_Name = ""
    for fasta in os.listdir(input_Folder):
        with open (f"{input_Folder}/{fasta}","r") as r:
            for line in r:
                if line.startswith(">"):
                    seq_Name = line.strip()
                    if seq_Name not in concatenatedFasta:
                        concatenatedFasta[seq_Name] = []
                else:
                    seq = line.strip()
                    concatenatedFasta[seq_Name].append(seq)
    with open ("concat.fasta","w") as f:
        for key, value in concatenatedFasta.items():
            f.write('%s\n' %(key))
            for i in value:
                f.write(i+"\n")

if __name__ == "__main__":
    getNamed_fastas('named_Fastas')
    fastas_Aligner('named_Fastas','alignedFastas')
    fastas_Concatenator('alignedFastas')

