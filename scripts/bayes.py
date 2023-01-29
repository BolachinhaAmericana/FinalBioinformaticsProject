import os

os.system("seqmagick convert --output-format nexus --alphabet dna concat.fasta concat.nex")
with open("concat.nex", "a") as f:
    f.write(""" 
begin mrbayes;
set autoclose=yes;
mcmc ngen=1000000 printfreq=500 samplefreq=500 diagnfreq=5000 nchains=4 savebrlens=yes;
mcmc;
smut filename=Myrun01;
end;
""")
os.system("mb concat.nex")
