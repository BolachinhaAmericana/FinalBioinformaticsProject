#!/usr/bin/env python3

import os
import subprocess

#This script will receive 2 inputs that will be one folder with all the fasta files and the second is the list of species names.


#concatenation of the fasta files
conca = ''
for file in os.listdir(pasta):
    with open(fastafile, 'r') as f:
        conca+= f.read()

print(conca)
