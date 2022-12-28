#!/usr/bin/env python3

import os
import subprocess

#This script will receive 2 inputs that will be one folder with all the fasta files and the second is the list of species names.


#concatenation of the fasta files
conca = ''
for file in os.listdir(pasta):
    with open(fastafile, 'r') as f:
        conca+= f.read()


#Alterar o nome do fasta concatenado
listaNomes= ["X1", "X2"]
contador = 0
with open("test.fasta", "r") as r:
    with open ("newtxt","w") as f:
        for line in r:
            if line.startswith(">"):
                f.write(">"+ listaNomes[contador] +"\n")
                contador +=1
            else:
                f.write(line)
        f.close()
