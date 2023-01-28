#!/usr/bin/env python3

import sys
import shutil
import os

ficheiro = sys.argv[1]


def nexusbasic(ficheiro):
	with open(f"{ficheiro}.nex", 'w+') as w:
		w.write(f"""#NEXUS\n
BEGIN DATA;
DIMENSIONS NTAX={len(dict)} NCHAR={len(dict[listaNomes[0]])};
FORMAT DATATYPE=DNA MISSING=N GAP=-;
MATRIX\n""")
		for key, value in dict.items():
			w.write(key +'  '+ value + '\n')
		w.write(f""" ;
END;\n
begin mrbayes;
set autoclose=yes;
mcmcp ngen=1000 printfreq=1000 samplefreq=100 diagnfreq=1000 nchains=4 savebrlens=yes filename=MyRun01;
mcmc;
sumt filename=MyRun01;
end;
""")
		return w


def nomeSeq(ficheiro):
	listaNomes =  []
	with open (f'{ficheiro}.fasta') as f:
		for lines in f:
			if '>' in lines:
				if len(lines) <100:
					linhas = lines.replace('\n', '')
					linha = linhas.replace('>','').replace(' ','')
					listaNomes.append(linha)
		return listaNomes


def seqConteudo(ficheiro):
	listaSeq =[]
	with open (f'{ficheiro}.fasta') as f:
		for lines in f:
			if '>' not in lines:
				linha = lines.replace('\n','')
				listaSeq.append(linha)
		return listaSeq

def obterSeqCorretas(listaSeq):
	seq_corretas =[]
	curr_seq = ''
	for i in range(len(listaSeq)):
		curr_seq = curr_seq + listaSeq[i]
		if len(listaSeq[i]) <= 59:
			seq_corretas.append(curr_seq)
			curr_seq = ''
			continue
	return seq_corretas

def obterDict(seq_corretas, Nomes):
	dict= {}
	dict = dict.fromkeys(Nomes)
	x = 0
	while x< len(dict):
		for l in listaNomes:
			dict[l] = seq_corretas[x]
			x = x+1
	return dict


def nexusToStdout(ficheiro):
	with open(f"{ficheiro}.nex", 'r') as r:
		shutil.copyfileobj(r, sys.stdout)

	return r

if __name__ ==  '__main__':
	listaNomes = nomeSeq(ficheiro)
	listaSeq = seqConteudo(ficheiro)
	seq_corretas = obterSeqCorretas(listaSeq)
	dict = obterDict(seq_corretas,listaNomes)
	nexusbasic(ficheiro)
	nexusToStdout(ficheiro)
