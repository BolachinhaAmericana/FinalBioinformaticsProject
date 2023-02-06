import os
from Bio import Entrez
from scripts.getSpeciesGeneList import env_setup
import shutil

def test_env_setup():
    output_dir = 'TestSuit/GeneLists'
    env_setup(output_dir)
    assert os.path.exists(output_dir)
    shutil.rmtree(output_dir)
