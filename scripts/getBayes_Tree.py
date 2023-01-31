#!/usr/bin/env python3
import subprocess
subprocess.run("seqmagick convert --output-format nexus --alphabet dna concat.fasta concat.nex", shell = True)
with open("concat.nex", "a") as f:
    f.write(""" 
begin mrbayes;
set autoclose=yes;
mcmc ngen=100000 printfreq=500 samplefreq=500 diagnfreq=5000 nchains=4 savebrlens=yes filename=MyBayes;
mcmc;
sump;
sumt conformat=simple contype=halfcompat;
end;
""")

subprocess.run("mb concat.nex", shell = True)
