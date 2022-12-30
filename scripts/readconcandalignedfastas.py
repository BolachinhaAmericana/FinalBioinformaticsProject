#!/usr/bin/env python3

import os
import subprocess

#This script will receive 2 inputs that will be one folder with all the fasta files and the second is the list of species names.

#alinhar fastas
def alignFastas(inputFolder,outputFolder)
    folder= inputFolder
    outfolder = outputFolder

    fasta_files = [f for f in os.listdir(folder) if f.endswith('.fasta')]


    for fasta in fasta_files:
        f = open(f"{outfolder}/_aligned{fasta}","w")
        subprocess.run(['mafft', folder + '/' + fasta], stdout = f)
    return outfolder

#concatenate all the fastas
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

    with open ("concaFasta.fasta","w") as f:
        for key, value in concatenatedFasta.items():
            f.write('%s\n' %(key))
            for i in value:
                f.write(i+"\n")


def getlistofNames(namesFile):
    listaNomes= []
    with open (f"{namesfile}", "r") as r:
        for line in r:
            listaNomes.append(line.strip())
        return listaNomes

def fastaFinal(listNames,concatenatedFasta):
    contador = 0
    with open(f"{concatenatedFasta}", "r") as r:
        with open ("concatenatedFinal.fasta","w") as f:
            for line in r:
                if line.startswith(">"):
                    f.write(">"+ listaNomes[contador] +"\n")
                    contador +=1
                else:
                    f.write(line)
            f.close()
