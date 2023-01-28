#!/usr/bin/env python3

import subprocess


def getBestModel(fasta):

    subprocess.run(f'modeltest-ng -d nt -i {fasta} -o model -p 4 -f e -h i -s 11 > model.txt', shell = True)

    with open("model.txt") as file:
        content = file.readlines()
        lastLine = content[-5]
    return lastLine

def raxMLExecutor(fasta, lastLine):
    if "+I" in lastLine:
        subprocess.run("raxmlHPC -s concat.fasta -n tree.tre -m GTRCAT -p 256789 -x 256789 -# autoFC", shell = True)
    else:
        subprocess.run("raxmlHPC -s concat.fasta -n tree.tre -m GTRCAT -p 256789 -x 256789 -# autoFC", shell = True)

if __name__ == "__main__":
    lastLine = ""
    getBestModel("concat.fasta")
    raxMLExecutor("concat.fasta", lastLine)
