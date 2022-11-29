#Doker commands:
#docker build -t imagename .
#docker run --name "conteiner name" -v $(pwd):/lab -it imagename

#Snakemake commands: 
#pip install snakemake==5.26.1
#python3 -m pip install snakemake
#term="Lagomorpha" rank="kingdom"  snakemake --cores all
#taxonkit list --ids 33090 -nr --indent "    "

import os
term = os.environ.get("term")
rank = os.environ.get("rank")


rule getTaxonIdRankTerm:
    params:
        cientificName = term,  
        rank = rank
    output:
        "taxonId.fasta"
    shell:
        "python3 scripts/getTaxonId.py {params.cientificName} {params.rank} > {output}"


rule getALLRankOrganism:
    output:
        "output.fasta"
    shell:
        'taxonkit list --ids $(cat taxonId.fasta) -nr --indent "    " > {output}'