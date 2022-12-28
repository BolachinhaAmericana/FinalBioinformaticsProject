#!/usr/bin/env python3

import subprocess

#-d input file
#-S Best Model
#-tr number of cores
#-p Criterion model
#-s number of substituion 11 to use all models defined by jmodeltest documentation
#-v model averaging and parameter importance
#-t for choose what the base tree. Maximum likelihoodL in this example
#fasta = ficheiro fasta concatenado recebido recebido

subprocess.run(f'jmodeltest -d {fasta} -tr 4  -p  -AIC -s 11 -i -v -t ML -o teste2', shell = True)

with open("teste2") as file:
    for line in file:
        pass
    last_line = line

if "+i" in last_line:
    print("Executar comando raxML com derivação +I")
else:
    print("Executar comando raxML sem derivação +I")


subprocess.run("rm teste2", shell = True)
