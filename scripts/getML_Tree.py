#!/usr/bin/env python3

import subprocess


def getBestModel(fasta):

    subprocess.run(f'modeltest-ng -d nt -i {fasta} -o model -p 4 -f e -h i -s 11 > model.txt', shell = True)

    with open("model.txt") as file:
        content = file.readlines()
        lastLine = content[-5]
    return lastLine

def raxMLExecutor(lastLine):
    if "+I" in lastLine:
        subprocess.run("raxmlHPC -m GTRCATI -p 1234567 -b 1234567 -# autoFC -s concat.fasta -n nhk", shell = True)
    else:
        subprocess.run("raxmlHPC -m GTRCAT -p 1234567 -b 1234567 -# autoFC -s concat.fasta -n nhk", shell = True)
    subprocess.run("rm -r model*", shell = True)
    subprocess.run("RAxML_info.nhk", shell = True) 

if __name__ == "__main__":
    lastLine = ""
    getBestModel("concat.fasta")
    raxMLExecutor(lastLine)
