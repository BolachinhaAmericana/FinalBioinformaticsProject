#!/usr/bin/env python3
'''
This file will output the taxon rank name of a input species.
'''

import sys
from Bio import Entrez

def get_user_arguments(): ######
    '''
    What for:
        Takes 2 arguments as std.input.
        These should be a scientific name of a species and a taxonomy tank

    Arguments:
        None

    Vars:
        scientific_name: Search term. Should be a species scientific name following the example: "Canis lupus".
        taxonomy_rank: Taxonomy Rank relation. If you dont know what to do here enter an item of this list:
            ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    Returns:
        scientific_name and taxonomy_rank
    '''
    try:
        scientific_name = sys.argv[1]
        taxonomy_rank = sys.argv[2]
        return scientific_name, taxonomy_rank
    except IndexError:
        sys.exit('''
                This program takes 2 arguments.
                1 - Species Scientific Name
                2 - Taxonomy rank (One item in the following list: 
                ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
''')

def get_taxon_name(scientific_name: str,taxonomy_rank: str): ######
    '''
    What for:
        Returns the scientific name of a taxon based on species scientific name and taxonomy rank.

    Args:
        scientific_name
        taxonomy_rank: Taxonomic rank of the taxon to search for. Must be one of the following:
            ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    
    Vars:
        search_handle: search query
        search_record: result of search query
        taxon_id: taxon id from our argument species
        fetch_handle: fetch query
        fetch_record: result of fetch query

    Returns:
        taxon rank with from the assigned species.
        Example: species= homo sapiens taxonomy_rank= order will return primates.
    '''
    Entrez.email = "no_warnings@email.com"
    search_handle = Entrez.esearch(db="taxonomy", term=scientific_name)
    search_record = Entrez.read(search_handle)

    taxon_id = search_record["IdList"][0]

    fetch_handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
    fetch_record = Entrez.read(fetch_handle)

    for taxon in fetch_record[0]["LineageEx"]:
        if taxon["Rank"] == taxonomy_rank:
            return taxon["ScientificName"]

def error_handle_taxon_name(scientific_name,taxonomy_rank): ######
    '''
    What for:
        Uses get_taxon_name. With this it will verify if the scientific name is correct.

    Vars:
        specie (str): Scientific name of the species to search for.
        rank (str): Taxonomic rank of the taxon to search for. Must be one of the following:
            ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            
    Returns:
        taxon
    '''
    try:
        taxon= get_taxon_name(scientific_name,taxonomy_rank)
        return taxon
    except IndexError:
        sys.exit("error: Bad Scientific Name")
        return None

if __name__ == "__main__":
    species,rank = get_user_arguments()
    taxon_name = error_handle_taxon_name(species, rank)
    if taxon_name:
        print(taxon_name)
    else:
        print(f"No taxon with rank {rank} was found for species {species}")
