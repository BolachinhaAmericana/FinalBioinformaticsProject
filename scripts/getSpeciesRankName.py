#!/usr/bin/env python3

#Made by Silva feat. Valente
from Bio import Entrez
import sys


def getUserArgs(specie,rank):
    '''
    Takes 2 arguments as std.input
    
    Arguments:
         None
    Vars -
        term: Search term. Should be a species scientific name following the example: "Canis lupus". Note that the 1st letter is the only capital.
        taxonRank: Taxonomic Rank relation. If you dont know what to put here just enter an item of the following list: 
            ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    Returns:
        term and taxonRank
    '''
    specie = sys.argv[1]
    rank = sys.argv[2] 
    return specie, rank

def get_taxon_name(specie,rank):
    '''
    Returns the scientific name of a taxon based on a search term and taxonomic rank.
    
    Vars:
        specie (str): Scientific name of the species to search for.
        rank (str): Taxonomic rank of the taxon to search for. Must be one of the following:
            ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    Returns:

        str: Scientific name of the taxon with the specified rank.
        None: If no taxon with the specified rank is found.
    '''
    try:
        Entrez.email = "no_warnings@email.com"
        handle = Entrez.esearch(db="taxonomy", term=specie)
        record = Entrez.read(handle)
        taxon_id = record["IdList"][0]
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        record = Entrez.read(handle)
        for taxon in record[0]["LineageEx"]:
            if taxon["Rank"] == rank:
                return taxon["ScientificName"]
    except IndexError:
        sys.exit("error: Bad ScientificName")
    return None

if __name__ == "__main__":
    specie,rank = getUserArgs(specie='',rank='')
    taxon_name = get_taxon_name(specie, rank)
    if taxon_name:
        print(taxon_name) 
    else:
        print("No taxon with rank {} was found for species {}".format(rank, specie))

