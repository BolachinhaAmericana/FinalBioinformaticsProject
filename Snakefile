#Doker commands:
#docker build -t imagename .
#docker run --name "conteiner name" -v $(pwd):/lab -it imagename

#Snakemake commands: 
#pip install snakemake==5.26.1
#python3 -m pip install snakemake
#dataBase="Nucleotide" term="animal name for example" snakemake --cores all

import os
dataBase = os.environ.get("dataBase")
term = os.environ.get("term")


rule getSecBdTerm:
    params:
        bd = dataBase,
        name = term    
    output:
        "output.fasta"
    shell:
        "python3 scripts/getSecBdTerm.py {params.bd} {params.name} > {output}"
   


