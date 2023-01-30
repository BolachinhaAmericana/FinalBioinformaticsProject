#!/usr/bin/env python3
import os

os.system("seqmagick convert --output-format nexus --alphabet dna concat.fasta concat.nex")
with open("concat.nex", "a") as f:
    f.write(""" 
begin mrbayes;
set autoclose=yes;
mcmc ngen=1000 printfreq=500 samplefreq=500 diagnfreq=5000 nchains=4 savebrlens=yes filename=MyRun01;
mcmc;
sump;
sumt;
output format=newick;
end;
""")

os.system("mb concat.nex")
os.system("mv MyRun01.con.tre MyBayes.tre")
os.system("rm -r MyRun*")
