#!/usr/bin/env python3


########################################### Test_Version ################################################################

import os
import pytest
from scripts.getSpeciesNamesList import downloadScientificNamesList
import scripts.essentials


def test_scientific_names_list_downloader():
    # test the function with a known input
    rank_name = 'primates'
    downloadScientificNamesList(rank_name)
    # check if the output file exists
    assert os.path.exists('ScientificNames_list.txt')
    # check if the file is not empty
    assert os.stat('ScientificNames_list.txt').st_size != 0

def test_scientific_names_list_downloader_without_input():
    # test the function with no input
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        downloadScientificNamesList()
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 'No input was give. This program takes a species taxonomy level (ex: primates).'
