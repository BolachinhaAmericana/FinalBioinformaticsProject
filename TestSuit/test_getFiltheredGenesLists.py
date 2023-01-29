import os
import pytest
import shutil
import sys

from scripts.getFiltheredGenesLists import get_user_arguments, proximity_tester


def test_get_user_arguments():
    with pytest.raises(SystemExit):
        get_user_arguments()
    sys.argv = ['script_name', 'search_term', '50', '70']
    assert get_user_arguments() == ('search_term', 50, 70)

def test_proximity_tester():
    target_list = {'gene1', 'gene2', 'gene3'}
    species_genes_list = {'gene1', 'gene4', 'gene5'}
    proximity_value = 50
    assert proximity_tester(species_genes_list, target_list, proximity_value) == False
    proximity_value = 25
    assert proximity_tester(species_genes_list, target_list, proximity_value) == True
