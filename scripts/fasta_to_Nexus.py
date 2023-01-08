#!/usr/bin/env python3

import sys
import shutil
import os

file = sys.argv[1]
ngen = sys.argv[2]
outgroup = sys.argv[3]


def nexusbasic(file):
	with open(f"{file}.nex", 'w+') as w:
		w.write(f"""#NEXUS\n
BEGIN DATA;
DIMENSIONS NTAX={len(dict)} NCHAR={len(dict[nameList[0]])};
FORMAT DATATYPE=DNA MISSING=N GAP=-;
MATRIX\n""")
		for key, value in dict.items():
			w.write(' '+ key +'  '+ value + '\n')
		w.write(f""" ;
END;\n
begin mrbayes;
set autoclose=yes;
outgroup {outgroup};
mcmcp ngen={ngen} printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename=MyRun01;
mcmc;
sumt filename=MyRun01;
end;
""")
		return w


def nameSeq(file):

	'''
    What for
	Get every name of gene seqs

    Arguments: 
        Fasta file

    Vars -
        nameList: list of names

    Returns:
        nameList it self
    '''

	nameList =  []
	with open (f'{file}.fasta') as f:
		for lines in f:
			if '>' in lines:
				if len(lines) <100:
					liness = lines.replace('\n', '')
					line = liness.replace('>','').replace(' ','')
					nameList.append(line)
		return nameList


def seqContent(file):

	'''
	What for
	Get the genes seqs it self

    Arguments: 
        Fasta file

    Vars -
        listSeq: all gene seqs

    Returns:
        listSeq it self

    '''

	seqList =[]
	with open (f'{file}.fasta') as f:
		for lines in f:
			if '>' not in lines:
				line = lines.replace('\n','')
				seqList.append(line)
		return seqList

def getCorrectSeq(seqList):

	'''
    What for
	Get only the wanted seqs

    Arguments: 
        seqList

    Vars -
        correct_seq: list with the wanted seqs

    Returns:
        correct_seq
    '''

	correct_seq =[]
	curr_seq = ''
	for i in range(len(seqList)):
		curr_seq = curr_seq + seqList[i]
		if len(seqList[i]) <= 59:
			correct_seq.append(curr_seq)
			curr_seq = ''
			continue
	return correct_seq

def getDict(correct_seq, Names):

	'''
    What for:
	Create a dictionary object where the keys are the names in a list and the values are the corresponding elements in another list. 
	
    Arguments: 
    	correct_seq: A list of elements to be used as values in the dictionary.
    	Names: A list of elements to be used as keys in the dictionary.

    Vars:
        dict: Dictionary that will store the keys and values.

    Returns:
        dict: Dictionary where the keys are the names in the 'Names' list and the values are the corresponding elements in the 'correct_seq' list.

    '''

	dict= {}
	dict = dict.fromkeys(Names)
	x = 0
	while x< len(dict):
		for l in nameList:
			dict[l] = correct_seq[x]
			x = x+1
	return dict

def mrBayes(outgroup):

	'''
    What for: 
	Searches if the outgroup correspond to a organism that was stored in nameList

    Arguments: 
        outgroup

    Vars -
        none

    Returns:
        outgroup 

    '''

	if outgroup not in nameList:
		print('Organism not found')
	else:
		return outgroup

def nexusToStdout(file):

	'''
    What for: 
	Get nexus file to standard output

    Arguments: 
        nexus file 

    Vars -
        none 

    Returns:
        output

    '''

	with open(f"{file}.nex", 'r') as r:
		shutil.copyfileobj(r, sys.stdout)

	return r

def deleteNexus(file):

	'''
	What for:
	Remove old nexus files

    Arguments: 
        Nexus file

    Vars -
        none

    Returns:
        none
    '''

	os.remove(f'{file}.nex')
	return

if __name__ ==  '__main__':
	nameList = nameSeq(file)
	seqList = seqContent(file)
	correct_seq = getCorrectSeq(seqList)
	dict = getDict(correct_seq,nameList)
	mrBayes(outgroup)
	nexusbasic(file)
	nexusToStdout(file)
	deleteNexus(file)
