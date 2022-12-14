import essentials as e 

from Bio import Entrez
import taxoniq
import os

def getTaxonId(eSearchResult): # 2.1
    '''
    Will get the Taxon Id of our input search term

    Arguments: 
        Takes eSearch result as argument.
    Vars - 
        taxonId: filthering string
    Returns:
        Filthered Taxon Id
    '''
    taxonId = str(eSearchResult["IdList"])
    taxonId = taxonId.replace("[\'",'')
    taxonId = taxonId.replace("\']",'')
    return taxonId

def getRankName(taxonId, taxonRank): #2.2
    '''
    Will get the specific name of the rank according to species and overall rank. Example: Human, order = Primates

    Arguments:
        Takes target taxon id and input taxonomy rank as arguments.
    Vars -
        t: loads taxon on var using Taxoniq API
        rankList: List of lineage of loaded Taxon
        relation: Compares the lineage with taxon ranks
        tax: taxon id of the related rank
        taxRankName: conversion of taxon id to scientific name using Taxoniq API
    Returns:
        taxRankName
    '''
    t = taxoniq.Taxon(taxonId)
    rankList = t.ranked_lineage
    relation = [(t.rank.name, t) for t in rankList]
    for i in relation:
        if taxonRank in str(i):
            tax = i[1]
            taxRankName = tax.scientific_name
            break
    return taxRankName
#1- Done- Taxonomy Rank Name (Primates)

def downloadScientificNamesList(RankName: str):
    '''
    Downloads ScientificNames_list.txt, with all species from said taxonomy rank name
    Args:
        RankName: Name of the taxonomy rank according to said specie. Example: Human order = Primates
    Vars-
        webEnv and queryKey: Web Environment and Query Key
        filther: filthering instance from efetch query.
    Returns: 
        None, but will create a ScientificNames_list.txt file with the scientific name of all species from RankName Argument.
    '''
    webEnv, queryKey= e.eSearch('Taxonomy', RankName+'[subtree]')["WebEnv"], e.eSearch('Taxonomy', RankName+'[subtree]')["QueryKey"]
    e.printOutputToFile('./eFetch.dump', 'w', e.eFetch('Taxonomy', webEnv, queryKey)) #get eFetch
    filther= e.readFile('./eFetch.dump').replace('<ScientificName>', '\n')
    e.printOutputToFile('./filthering.dump', 'w', filther)
    os.system("grep 'species' filthering.dump > list.dump")
    os.system("grep -v 'subspecies' list.dump > loading.dump")
    filther= e.readFile('./loading.dump').replace('</ScientificName>', '\n')
    e.printOutputToFile('./finalFilther.dump', 'w', filther)
    os.system("grep -v 'species' ./finalFilther.dump > ./ScientificNames_list.txt")


    os.system('rm ./*.dump')
    return None
#2- Done- ScientificNames_list.txt
def downloadSpecieGeneList(specie: str, batch_size: int): #4.1
    '''
    Downloads one species gene list
    Arguments:
        specie: scientific name of a specie to get gene list
        batch_size: max= 10000. The number of genes that it will fetch at a time.
    Vars-
        eSearchResults: results of a eSearch on the Gene db, filthering all replaced or discontinued genes.
        count: number of total genes that will be fetched
        webenv: Web Environment
        query_key: Query Key
        fetch_handle: eFetch Query.
        data: read of the fetch_handle
    Returns: None, but will create a LoadingGeneList.txt file

    '''
    eSearchResults= e.eSearch('Gene', f'{specie} AND alive[prop]')
    count = int(eSearchResults["Count"])
    webenv = eSearchResults["WebEnv"]
    query_key = eSearchResults["QueryKey"]

    for start in range(0, count, batch_size):
        while True:
            try:
                print(f"Downloading record {start+1} to {count} - {specie}")
                # Can't use eFetch function because if has a lot of extra paramethers and might aswell just leave it be.
                fetch_handle = Entrez.efetch(db="Gene",rettype="gb",retmode="text",retstart=start,retmax=batch_size,webenv=webenv,query_key=query_key,idtype="acc",)
                data = fetch_handle.read()
                fetch_handle.close()
                e.printOutputToFile('loading.dump', 'w', data)
                # for some reason isn't accepting subprocess. So we used os instead ,that works fine.
                #subprocess.call('./bashScripts/getGeneListFilthered.sh')
                
                os.system("grep -v ':' loading.dump > loading1.dump")
                os.system("grep '\.' loading1.dump > loading2.dump")
                os.system("grep -v '\[' loading2.dump > loading3.dump")
                os.system("sed 's/^.*\ //' loading3.dump >> LoadingGeneList.txt")
                os.system('rm ./*.dump')
                break
            except:
                print('Error while downloading. Trying again...')
    return None #Outputs the file  LoadingGeneList.txt

def organizeSpecieGeneList(specie: str, pathToGeneList):
    '''
    Organizes the gene list and reacts incase list is empty.
    Arguments:
        takes the specie name as an argument
    Vars-
        None
    Returns:
        None, but will rename and move a file to GeneLists/
    '''
    try:
        os.rename(pathToGeneList, './GeneLists/{}_GeneList'.format(specie))
    except:
        e.printOutputToFile('./GeneLists/noGeneSpecies.txt', 'a', specie)
        print(f'The following species has zero genes: {specie}')
    return None
#3- Done- Genes Lists








term, taxonRank= e.getUserArgs()
RankName= getRankName(getTaxonId(e.eSearch('Taxonomy', term)), taxonRank)
#1- Done
downloadScientificNamesList(RankName)
#2- Done
speciesList= e.readFileToList('./ScientificNames_list.txt')
if os.path.exists('./GeneLists'):
    os.system('rm ./GeneLists/*')

else: os.system('mkdir ./GeneLists')
for specie in speciesList:
    downloadSpecieGeneList(specie, 7500)
    organizeSpecieGeneList(specie, './LoadingGeneList.txt')
#3- Done
