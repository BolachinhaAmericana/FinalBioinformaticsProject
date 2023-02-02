#!/usr/bin/env python3

import subprocess


def getBestModel(fasta):
    '''
    What for:
        Runs modeltest-ng on a fasta file to determine the best model of evolution.

    Arguments:
        fasta: Path to the fasta file to run modeltest-ng on.
    
    Vars:
        content: list that contains the lines of the file "Model.txt".
        lastLine: Variable that is used to store the last line of the file "Model.txt".

    Returns:
       lastLine

    '''
    subprocess.run(f'modeltest-ng -d nt -i {fasta} -o model -p 4 -f e -h i -s 11 > model.txt', shell = True)

    with open("model.txt") as file:
        content = file.readlines()
        lastLine = content[-5]
    return lastLine

def raxMLExecutor(lastLine):
    '''
    What for:
        This function is used to run RaxML on a given FASTA file.

    Arguments/Vars:
        fasta: Path of the FASTA file to be used as input for RaxML.
        lastLine: The last line of the file.
        
    Returns:
       None

    '''
    if "+I" in lastLine:
        subprocess.run("raxmlHPC -m GTRCATI -p 1234567 -x 1234567 -# autoFC -s concat.fasta -n nwk", shell = True)
    else:
        subprocess.run("raxmlHPC -m GTRCAT -p 1234567 -x 1234567 -# autoFC -s concat.fasta -n nwk", shell = True)
    #subprocess.run("rm -r model*", shell = True)
    #subprocess.run("RAxML_info.nwk", shell = True) 

if __name__ == "__main__":
    lastLine = ""
    getBestModel("concat.fasta")
    raxMLExecutor(lastLine)
