#!/usr/bin/env python3


########################################### Test_Version ################################################################
import sys
from scripts.getSpeciesRankName import get_user_arguments, get_taxon_name,error_handle_taxon_name

def test_get_user_arguments():
    sys.argv = ["getSpeciesRankName.py","Passer domesticus", "family"]
    scientific_name, taxonomy_rank = get_user_arguments(sys.argv[1], sys.argv[2])
    assert scientific_name == "Passer domesticus"
    assert taxonomy_rank == "family"

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