#!/usr/bin/env python3

import subprocess

#-d input file
#-S Best Model
#-tr number of cores
#-p Criterion model
#-s number of substituion 11 to use all models defined by jmodeltest documentation
#-v model averaging and parameter importance
#-t for choose what the base tree. Maximum likelihoodL in this example
#fasta = ficheiro fasta concatenado recebido

#subprocess.run(f'modeltest-ng -d nt, shell = True)

with open("help.txt") as file:
    content = file.readlines()
    lastLine = content[-4]
if "+I" in lastLine:
   print("Executar comando raxML com derivação +I")
else:
    print("Executar comando raxML sem derivação +I")
