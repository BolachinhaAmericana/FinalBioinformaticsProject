#!/usr/bin/env python3
from Bio import SeqIO
import subprocess
import shutil
import glob
import os
# This script will receive 2 inputs that will be one folder with all the fasta files and the second is the list of species names.

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


# Alinhar fastas
def alignFastas(inputFolder,outfolder):
    folder = inputFolder
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    fasta_Files = [f for f in os.listdir(folder) if f.endswith('.fasta')]
    for fasta in fasta_Files:
        f = open(f"{outfolder}/_aligned{fasta}","w")
        subprocess.run(['mafft', folder + '/' + fasta], stdout = f)


# Concatenation of aligned fasta files
def concatenateFastas(inputFolder):
    folder = inputFolder
    concatenatedFasta = {}
    seq = ""
    name = ""
    for fasta in os.listdir(folder):
        with open (f"{folder}/{fasta}","r") as r:
            for line in r:
                if line.startswith(">"):
                    name = line.strip()
                    if name not in concatenatedFasta:
                        concatenatedFasta[name] = []
                else:
                    seq = line.strip()
                    concatenatedFasta[name].append(seq)
    subprocess.run(f'rm -r {folder}', shell = True)
    with open ("concat.fasta","w") as f:
        for key, value in concatenatedFasta.items():
            f.write('%s\n' %(key))
            for i in value:
                f.write(i+"\n")

if __name__ == "__main__":
    getNamed_fastas('named_Fastas')
    alignFastas('named_Fastas','alignedFastas')
    concatenateFastas('alignedFastas')
