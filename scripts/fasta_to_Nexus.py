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
	nameList =  []
	with open (f'{file}.fasta') as f:
		for lines in f:
			if '>' in lines:
				if len(lines) <100:
					linhas = lines.replace('\n', '')
					linha = linhas.replace('>','').replace(' ','')
					nameList.append(linha)
		return nameList


def seqContent(file):
	listSeq =[]
	with open (f'{file}.fasta') as f:
		for lines in f:
			if '>' not in lines:
				linha = lines.replace('\n','')
				listSeq.append(linha)
		return listSeq

def getCorrectSeq(listaSeq):
	correct_seq =[]
	curr_seq = ''
	for i in range(len(listaSeq)):
		curr_seq = curr_seq + listaSeq[i]
		if len(listaSeq[i]) <= 59:
			correct_seq.append(curr_seq)
			curr_seq = ''
			continue
	return correct_seq

def obterDict(correct_seq, Nomes):
	dict= {}
	dict = dict.fromkeys(Nomes)
	x = 0
	while x< len(dict):
		for l in nameList:
			dict[l] = correct_seq[x]
			x = x+1
	return dict

def mrBayes(outgroup):
	if outgroup not in nameList:
		print('Organism not found')
	else:
		return outgroup

def nexusToStdout(file):
	with open(f"{file}.nex", 'r') as r:
		shutil.copyfileobj(r, sys.stdout)

	return r

def deleteNexus(file):
	os.remove(f'{file}.nex')
	return
if __name__ ==  '__main__':
	nameList = nameSeq(file)
	seqList = seqContent(file)
	correct_seq = getCorrectSeq(seqList)
	dict = obterDict(correct_seq,nameList)
	mrBayes(outgroup)
	nexusbasic(file)
	nexusToStdout(file)
	deleteNexus(file)
