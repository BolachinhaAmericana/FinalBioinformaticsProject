#!/usr/bin/env python3

import subprocess


def getBestModel(fasta):

    subprocess.run(f'modeltest-ng -d nt -i {fasta} -o Model -p 4 -f e -h i -s 11 > teste2.txt', shell = True)

    with open("Model.txt") as file:
        content = file.readlines()
        lastLine = content[-5]
    return lastLine

def raxMLExecutor(fasta, lastLine):
    if "+I" in lastLine:
        subprocess.run("raxmlHPC -s concat.fasta -n teste -m GTRCAT -p 256789 -x 256789 -# autoFC", shell = True)
    else:
        subprocess.run("raxmlHPC -s concat.fasta -n teste -m GTRCAT -p 256789 -x 256789 -# autoFC", shell = True)

if __name__ == "__main__":
    lastLine = ""
    getBestModel("concat.fasta")
    executeRaxML("concat.fasta", lastLine)
