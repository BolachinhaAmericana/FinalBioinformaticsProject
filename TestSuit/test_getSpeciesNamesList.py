import os
import pytest
from scripts.getSpeciesNamesList import scientific_names_list_downloader

def test_scientific_names_list_downloader():
    # test the function with a known input
    rank_name = 'Passeridae'
    scientific_names_list_downloader(rank_name)
    # check if the output file exists
    assert os.path.exists('ScientificNames_list.txt')
    # check if the file is not empty
    assert os.stat('ScientificNames_list.txt').st_size != 0

