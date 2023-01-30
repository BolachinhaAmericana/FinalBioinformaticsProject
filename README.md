# Automatic Phylogenetic Tree

In this project, we aim to automate the creation and visualization of a phylogenetic tree based on a species input by the user. Our goal is to make the process of creating such trees more efficient and accessible, allowing for for deeper insights into evolutionary relationships within the animal kingdom.

## Getting Started

In order to run this project you will need to have Docker and fetch our image that can be found [here]

### Dependencies
The following dependencies are required to run this project:

> Docker v20.10.22\
snakemake v()
# Quando tiver acabado meter a versao do snakemake

### Installing

## Running the tests
# Need to do
Instructions on how to run any automated tests for the project.

## Usage
This program works using a set of files working with input and outputs 
# Need to do
Instructions on how to use the project, including any necessary configuration or settings.

## Authors
* **Duarte Valente** - *genes gathering and filtering* - [GitHub](https://github.com/BolachinhaAmericana) - [GitLab](https://gitlab.com/BolachinhaAmericana)
* **Gonçalo Alves** - *code reviwer and report manager* - - [GitHub](https://github.com/GonaloAlves) - [GitLab](https://gitlab.com/alvesgoncas2014)
* **Guilherme Silva** - *general manager* - [GitHub](https://github.com/GuilhermeVCCdaSilva) - [GitLab](https://gitlab.com/guilherme.vcc.silva)
* **Marine Fournier** - *overall helper and testing* - [GitHub](https://github.com/MarineF22) - [GitLab](https://gitlab.com/marine.fournier2002)
* **Matilde Machado** - *mr bayes tree and report writing* - [GitLab](https://gitlab.com/matildemachado)
* **Rodrigo Pinto** - *fasta management and maximum likelyhood tree* - - [GitHub](https://github.com/Sepay) - [GitLab](https://gitlab.com/Sepay1)

## Acknowledgments
- ChatCPT
We used external help from ChatGPT, an AI language model developed by OpenAI to improve our code and help with tests.

## How the programs work individually

Incase it is desired to run the programs individually and outside of the provided docker image, run the scripts in the following order:
1. getSpeciesRankName.py\
Takes a specie name and a taxonomy rank as arguments, returns the name of the especific taxon rank name. If this is confusing here is an example of how the code would work:
```bash
python3 getSpeciesRankName.py <specie_name> <taxonomy_rank>
python3 getSpeciesRankName.py 'Homo sapiens' 'order'
>> Primates
```
2. getSpeciesNamesList.py
Takes the output of the previous file, ```Primates``` as in the example and queries all species that share this taxonomy rank. Outputs a list of these species to ```ScientificNames_list.txt``` file.
3. getSpeciesGeneList.py
Runs the ```ScientificNames_list.txt``` file to check all species and downloads all gene names from each specie to ```/GeneLists/{specie}_GeneList```
4. getSatisfiedList.py
Filthers unecessary gene lists using two methods.
    - species filthering
will remove all species gene lists that don't make the cut of >input< percentage. Essensialy it compares the number of similar genes regarding the original species and if it shares a number above the (%), gets approved and clears the '1st stage'.
    - gene filthering
    In order to avoid the count of genes exclusive to the input species, we filther the genes that don't appear on a input (%) of other species that already passed the 'first stage'. Those who pass this 'second stage', get labeled as approved and are sent to ```FiltheredGeneNames_list.txt```