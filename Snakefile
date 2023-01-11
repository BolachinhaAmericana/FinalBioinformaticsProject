#Made by Silva

#Doker commands:
#docker build -t imagename .
#docker run --name "conteiner name" -v $(pwd):/lab -it imagename

#Snakemake commands: 
#pip install snakemake==5.26.1
#python3 -m pip install snakemake
#term="'Passer domesticus'" rank="family" proximity="30" similarity="30"  snakemake --cores all

import os
term = os.environ.get("term")
rank = os.environ.get("rank")
proximity = os.environ.get("proximity")
similarity = os.environ.get("similarity")

rule all:
    input:
        "RankName.txt",
        "ScientificNames_list.txt",
        directory("GeneLists"),
        "FiltredScientificNames_list.txt",
        "FiltredGeneNames_list.txt",
        directory("Squences_Fasta"),
        directory("named_Fastas"),
        "concat.fasta"      
           


rule getSpeciesRankName:
    params:
        cientificName = term,
        rank = rank
    output:
        "RankName.txt"
    shell:
        """
        python3 scripts/getSpeciesRankName.py {params.cientificName} {params.rank} > {output}
        """


rule getSpeciesNamesList:
    input:
        rules.getSpeciesRankName.output
    output:
        "ScientificNames_list.txt"
    shell:
        """
        python3 scripts/getSpeciesNamesList.py $(cat {input}) > {output}
        """
        

rule getSpeciesGeneList:
    input:
        rules.getSpeciesNamesList.output
    output:
        directory("GeneLists")
    shell:
        """
        python3 scripts/getSpeciesGeneList.py  
        """

rule getSatisfiedList:
    input:
        rules.getSpeciesGeneList.output
    params:
        cientificName = term,
        proximity = proximity,
        similarity = similarity

    output:
        "FiltredScientificNames_list.txt",
        "FiltredGeneNames_list.txt"
        
    shell:
        """
        python3 scripts/getSatisfiedList.py {params.cientificName} {params.proximity} {params.similarity}
        """ 

rule getSecToAlign:
    input:
        rules.getSatisfiedList.output
    output:
        directory("Squences_Fasta")
    shell:
        """
        python3 scripts/getSecToAlign.py
        """         

rule getConcatAlignNamedFasta:
    input:
        rules.getSecToAlign.output
    output:
        directory("named_Fastas"),
        "concat.fasta" 
    shell:
        """
        python3 scripts/getConcatAlignNamedFasta.py
        """  