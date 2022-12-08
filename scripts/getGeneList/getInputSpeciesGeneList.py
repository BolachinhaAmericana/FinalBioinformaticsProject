from Bio import Entrez
import sys
import os

f = open('eFetch.xml', 'r')
file = f.read()

#LIMITADO: consegue apernas ler os priemeiros 10k+- genes, porque tem limite de caracteres na mesma linha (.xml) 

first = file.replace('class="desc">', '\n') #Working!
second = first.replace('</p>', '\n')
with open('test.xml', 'w') as sys.stdout:
    print(second)

os.system('grep -v "div" test.xml > geneList.txt')
os.system('rm test.xml')

