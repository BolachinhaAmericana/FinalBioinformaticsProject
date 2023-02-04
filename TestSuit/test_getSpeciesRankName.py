#!/usr/bin/env python3


########################################### Test_Version ################################################################
import sys
from scripts.getSpeciesRankName import  get_taxon_name,error_handle_taxon_name

def test_get_taxon_name():
    scientific_name = "Passer Domesticus"
    taxonomy_rank = "family"
    taxon_name = get_taxon_name(scientific_name, taxonomy_rank)
    assert taxon_name == "Passeridae"

def test_error_handle_taxon_name():
    scientific_name = "Passer Domesticus"
    taxonomy_rank = "family"
    taxon_name = error_handle_taxon_name(scientific_name, taxonomy_rank)
    assert taxon_name == "Passeridae"