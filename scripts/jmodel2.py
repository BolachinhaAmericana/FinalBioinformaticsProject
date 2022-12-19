#!/usr/bin/env python3

import subprocess

#-d input file
#-S Best Model
#-tr number of cores
#-p Criterion model
#-s number of substituion 11 to use all models defined by jmodeltest documentation
#-v model averaging and parameter importance
#-t for choose what the base tree. Maximum likelihoodL in this example

subprocess.run('jmodeltest -d AK_s7rp_aln.fasta -S BEST -tr 4  -p  -AIC -s 11 -i -v -t ML -o teste2', shell = True)


