#!/usr/bin/env python3
'''mr bayes'''
import subprocess
subprocess.run("seqmagick convert --output-format nexus --alphabet dna concat.fasta concat.nex", shell = True, check = False)

with open("concat.nex", "a", encoding = 'utf8') as doc:
    doc.write("""
begin mrbayes;
set autoclose=yes;
mcmc ngen=100000 printfreq=500 samplefreq=500 diagnfreq=5000 nchains=4 savebrlens=yes filename=MyBayes;
mcmc;
sump;
sumt conformat=simple contype=halfcompat;
end;
""")
subprocess.run("mb concat.nex", shell = True, check = False)
