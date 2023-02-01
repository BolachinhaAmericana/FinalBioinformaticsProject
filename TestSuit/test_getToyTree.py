import sys
from scripts.get_toy_tree import get_user_input, get_user_args

def test_get_user_input():
    sys.argv = ['scriptname', 'tree1', 'tree2']
    MLTree, MrBayes = get_user_input(MLTree='', MrBayes='')
    assert MLTree == 'tree1'
    assert MrBayes == 'tree2'

def test_get_user_args():
    sys.argv = ['scriptname', 'tree1', 'tree2', 'species name']
    specieName = get_user_args(specieName='')
    assert specieName == 'species_name'

