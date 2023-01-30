#Made by Silva

#Doker commands:
#docker build -t imagename .
#docker run --name "conteiner name" -v $(pwd):/lab -it imagename

#Snakemake commands: 
#pip install snakemake==5.26.1
#python3 -m pip install snakemake
#term="'Passer domesticus'" rank="family" proximity="30" similarity="30"  snakemake --cores all

import os
import shutil

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
        "named_Fastas",
        "concat.fasta",
        directory("output_dir"),
        "RAxML_bootstrap.nhk", 
        "MyBayes.tre" 

rule copy_outputs:
    input:
        "RankName.txt",
        "ScientificNames_list.txt",
        "GeneLists",
        "FiltredScientificNames_list.txt",
        "FiltredGeneNames_list.txt",
        "Squences_Fasta",
        "named_Fastas",
        "concat.fasta",
        "RAxML_bootstrap.nhk", 
        "MyBayes.tre"
    output:
        directory("output_dir")
    run:
        if not os.path.exists("output_dir"):
            os.mkdir("output_dir")
        for item in input:
            if os.path.isdir(item):
                shutil.copytree(item, "output_dir"+"/"+os.path.basename(item))
            else: 
                shutil.copy(item, "output_dir")                
         
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

rule getMl_Tree:
    input:
        rules.getConcatAlignNamedFasta.output
    output:
        "RAxML_bootstrap.nhk"
    shell:
        """
        python3 scripts/getML_Tree.py 
        """

rule getBayes_Tree:
    input:
        rules.getMl_Tree.output
    output:
        "MyBayes.tre"
    shell:
        """
        python3 scripts/getBayes_Tree.py
        """

rule get_toy_tree:
    input:
        rules.getConcatAlignNamedFasta.output
    output:
        "RAxML_bootstrap.nhk"
        "MyBayes.tre"
    shell:
        """
        python3 scripts/get_toy_tree.py {params.cientificName}
        """

rule clean:
    input:
        "RankName.txt",
        "ScientificNames_list.txt",
        directory("GeneLists"),
        "FiltredScientificNames_list.txt",
        "FiltredGeneNames_list.txt",
        directory("Squences_Fasta"),
        "named_Fastas",
        "concat.fasta",
        "RAxML_bootstrap.nhk",
        "MyBayes.tre"
    shell:
        "rm -rf {input}"  

