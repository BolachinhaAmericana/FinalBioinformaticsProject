# Automagic phylogenies

In this project, we aim to automate the creation and visualization of a phylogeny tree based on a species input by the user. Our goal is to make the process of creating such trees more efficient and accessible, allowing for for deeper insights into evolutionary relationships within the animal kingdom.

## Dependencies
The following dependencies are required to run this project:
> Docker v20.10.22
- All other dependencies are automatically installed on our snake-make image. Therefore the only program the user **HAS** to have is docker.

## Running the tests
The tests are run every time there is a commit on the master branch of our repository.

## Usage
This program works using a set of files working with input and outputs. However you can **run** the program following these simple steps:

```bash
$ git clone https://gitlab.com/LabBinf2022/team-b/auto_magic_phylogenetic.git
# This line downloads the magic_phylogenetic directory with our project
$ cd auto_magic_phylogenetic
$ bash executer.sh -s <term> -r <rank> -p <proximity> -i <similarity>
```

### FAQ
#### Q: "Please provide all required arguments" Error, why is this happening?
##### Missing Brackets
You need to put all inputs between brackets ("") for the program to work
#### Q: How do I set the parameter correctly?
The variables that you need to pass to the program are the following:
##### Species (term)
Input the scientific name of the species that is intended for the study.
##### Taxonomy Rank (rank)
Input of the taxonomy level which complements the specie.\
In you unfamiliar with these levels, they can be:
- `kingdom`
- `phylum`
- `class`
- `order`
- `family`
- `genus`

Take in consideration that the recommended input on this parameter is `family` for the data to be neither too big nor too small for a decent and fast analysis.
##### Proximity Value
Percentage of proximity value. This parameter is better explained further down on this file!
##### Similarity Value
Percentage of similarity value. This parameter is better explained further down on this file!

## How the programs work individually

In-case it is desired to run the programs individually and outside of the provided docker image, run the scripts as following. `**Fair Warning**` keep in mind that you would have to install all additional programs that run on automatically on our snake-make image manually!!

1. **getSpeciesRankName.py**\
Takes a specie name and a taxonomy rank as arguments, returns the name of the specific taxon rank name.\
If this is confusing here is an example on how the code would work:

```bash
python3 getSpeciesRankName.py <specie_name> <taxonomy_rank>
python3 getSpeciesRankName.py 'Passer domesticus' 'order'
>> Passeriformes
# Another example would be that the Homo sapiens order is "Primates".
```

2. **getSpeciesNamesList.py**
Takes the output of the previous file, ```Passeriformes``` as in the example and queries all species that share this taxonomy rank. Outputs a list of these species to ```ScientificNames_list.txt``` file.

```bash
python3 getSpeciesNamesList.py <taxonomy_rank_name>
python3 getSpeciesNamesList.py 'Passeriformes'
$ ls
>> essentials.py            getSpeciesRankName.py    getSpeciesNamesList.py
   ScientificNames_list.txt
$ head -n 3 ScientificNames_list.txt 
>> 
Passer motitensis
Passer eminibey
Passer sp.
```

3. **getSpeciesGeneList.py**
Runs the ```ScientificNames_list.txt``` file to check all species and downloads all gene names from each specie to ```/GeneLists/{specie}_GeneList```

```bash
python3 getSpeciesGeneList.py <Takes 0 Outputs>
python3 getSpeciesGeneList.py
$ ls
>>  essentials.py    getSpeciesRankName.py    getSpeciesNamesList.py      
    GeneLists        ScientificNames_list.txt
$ cd GeneLists
$ ls
>> noGene   'Passer montanus_GeneList'      'Passer ammodendri_GeneList'
$ head -n 3 'Passer montanus_GeneList'
>> 
AZ368_gt21
ND2
ATP8
```

4. **getSatisfiedList.py**
Filters unnecessary gene lists using two methods.
    - species filtering (proximity)
will remove all species gene lists that don't make the cut of >input< percentage. Essentially it compares the number of similar genes regarding the original species and if it shares a number above the (%), gets approved and clears the '1st stage'.
    - gene filtering ()
    In order to avoid the count of genes exclusive to the input species, we filter the genes that don't appear on a input (%) of other species that already passed the 'first stage'. Those who pass this 'second stage', get labeled as approved and are sent to ```FiltheredGeneNames_list.txt```
    
```bash
python3 getSatisfiedList.py.py <input_specie> <proximity_value> <similarity_value>
python3 getSatisfiedList.py.py 'Passer domesticus' '10' '10'
$ ls
>>  essentials.py    getSpeciesRankName.py    getSpeciesNamesList.py      
    GeneLists        ScientificNames_list.txt
$ cd GeneLists
$ ls
>> noGene   'Passer montanus_GeneList'      'Passer ammodendri_GeneList'
$ head -n 3 'Passer montanus_GeneList'
>> 
AZ368_gt21
ND2
ATP8
```

5. **getSecToAlign.py** gathers genes from the filtered list of genes. Will take one by one and download all species outputs files: gene_name.fasta on the **Sequences_Fasta**

```bash
python3 getSecToAlign.py 
$ ls
>>  essentials.py    getSpeciesRankName.py    getSpeciesNamesList.py      
    GeneLists        ScientificNames_list.txt Squences_Fasta
$ cd Sequences_Fasta
$ ls
>> ATP6.fasta  COX1.fasta  COX3.fasta  ND1.fasta  ND3.fasta  ND4L.fasta  ND6.fasta
   ATP8.fasta  COX2.fasta  CYTB.fasta  ND2.fasta  ND4.fasta  ND5.fasta
# looking at one of the files
$ cat ATP6.fasta
>> 
>MN356394.1 Passer domesticus mitochondrion, complete genome
GTCCTTGTAGCTTAAAAAAAGCATGACACTGAAGATGTCAAGATGGCTGCCACACACACCCAAGGACAAA
...
>NC_029344.1 Passer ammodendri mitochondrion, complete genome
GTCCTGTAGCCTTATAAAAAGCATGACACTGAAGATGTCAAGATGGCTGCTACACGCACCCAAGGACAAA
...
```

6. **getConcatAlignNamedFasta.py** will align, rename the sequences on the files and concatenate the .fasta documents. Will output the named_Fastas directory with the .fasta with the updated names and the concat.fasta file with the concatenated fasta files.

```bash
python3 getConcatAlignNamedFasta.py
$ ls
>>  essentials.py    getSpeciesRankName.py    getSpeciesNamesList.py      
    GeneLists        ScientificNames_list.txt Squences_Fasta
    named_Fastas     concat.fasta
$ head -n 3 concat.fasta
>>
>Passer_ammodendri
gtcctgtagccttataaaaagcatgacactgaagatgtcaagatggctgctacacgcacc
caaggacaaaagacttagtcctaaccttactgttagtttttgctaggtatatacatgcaa
$ cd named_Fastas
$ ls
>> ATP6_updated.fasta  COX1_updated.fasta  COX3_updated.fasta  ND1_updated.fasta  ND3_updated.fasta  
head -n 3 ATP6_updated.fasta
>>
>Passer_ammodendri
GTCCTGTAGCCTTATAAAAAGCATGACACTGAAGATGTCAAGATGGCTGCTACACGCACC
CAAGGACAAAAGACTTAGTCCTAACCTTACTGTTAGTTTTTGCTAGGTATATACATGCAA
```

7. **getML_Tree.py** creates the maximum likelyhood tree using the concat.fasta file.

```bash
python3 getML_tree.py
>>/-Pyrgilauda_ruficollis
  |
  |   /-Onychostruthus_taczanowskii
  |--|
  |  |   /-Pyrgilauda_davidiana
  |   -|
--|      -Pyrgilauda_blanfordi
  |
  |      /-Montifringilla_adamsi
  |   /-|
  |  |  |   /-Montifringilla_nivalis
  |  |   -|
  |  |      -Montifringilla_henrici
   -|
     |      /-Passer_domesticus
     |   /-|
     |  |   -Passer_ammodendri
     |  |
      -|      /-Chloebia_gouldiae
        |   /-|
        |  |   -Padda_oryzivora
         -|
           |   /-Prunella_himalayana
            -|
              |   /-Prunella_rubeculoides
               -|
                 |   /-Prunella_strophiata
                  -|
                    |   /-Prunella_montanella
                     -|
                        -Prunella_fulvescens
```

8. **getBayes_tree.py** creates the mr bayes tree using concat.fasta file.
```bash
python3 getBayes_tree.py
>> /---------------------------------------------------------- Pyrgilauda_rufi~ (1)
   |                                                                               
   |                                   /---------------------- Prunella_rubecu~ (3)
   |                                   |                                           
   |                                   |              /------- Prunella_stroph~ (4)
   |                            /--100-+       /--100-+                            
   |                            |      |       |      \------- Prunella_monta~ (12)
   |                            |      \--100--+                                   
   |                     /--100-+              \-------------- Prunella_fulve~ (14)
   |                     |      |                                                  
   |                     |      \----------------------------- Prunella_himal~ (15)
   |              /--100-+                                                         
   |              |      |                            /------- Chloebia_gouldi~ (6)
   |              |      \-------------100------------+                            
   +      /--100--+                                   \------- Padda_oryzivora (10)
   |      |       |                                                                
   |      |       |                                   /------- Passer_ammoden~ (11)
   |      |       \----------------100----------------+                            
   |      |                                           \------- Passer_domesti~ (13)
   |--100-+                                                                        
   |      |                                           /------- Montifringilla_~ (8)
   |      |                                    /--100-+                            
   |      |                                    |      \------- Montifringilla_~ (9)
   |      \-----------------100----------------+                                   
   |                                           \-------------- Montifringilla~ (16)
   |                                                                               
   |                                           /-------------- Onychostruthus_~ (2)
   |                                           |                                   
   \--------------------100--------------------+      /------- Pyrgilauda_davi~ (5)
                                               \--100-+                            
                                                      \------- Pyrgilauda_blan~ (7)
```
9. **get_toy_tree.py** contains and defines information about the that are presented. Will take as input the tree files as well as the initial input species and create the final .pdf file containing the phylogeny trees from the analysis.
```bash
python3 get_toy_tree.py <ml_tree_file> <mb_tree_file> <input_species>
```

## Authors
* **Duarte Valente** - *genes gathering and filtering* - [GitHub](https://github.com/BolachinhaAmericana) - [GitLab](https://gitlab.com/BolachinhaAmericana)
* **Gonçalo Alves** - *code reviwer and report manager* - - [GitHub](https://github.com/GonaloAlves) - [GitLab](https://gitlab.com/alvesgoncas2014)
* **Guilherme Silva** - *project manager* - [GitHub](https://github.com/GuilhermeVCCdaSilva) - [GitLab](https://gitlab.com/guilherme.vcc.silva)
* **Marine Fournier** - *overall helper and testing* - [GitHub](https://github.com/MarineF22) - [GitLab](https://gitlab.com/marine.fournier2002)
* **Matilde Machado** - *mr bayes tree and report writing* - [GitLab](https://gitlab.com/matildemachado)
* **Rodrigo Pinto** - *fasta management and maximum likelyhood tree* - - [GitHub](https://github.com/Sepay) - [GitLab](https://gitlab.com/Sepay1)

## Acknowledgments
- ChatGPT-3
We used external help from ChatGPT, an AI language model developed by OpenAI to improve our code and help with tests. We also recurred to this tool to help with writing the report and this README!
