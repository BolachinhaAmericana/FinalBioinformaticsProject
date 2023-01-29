import os
import pytest
from scripts.getSecToAlign import createSquences_Fasta_Dir, getFastas

def test_createSquences_Fasta_Dir():
    directory = 'Squences_Fasta'
    pathToFastaDir = createSquences_Fasta_Dir(directory)
    assert os.path.exists(pathToFastaDir) == True