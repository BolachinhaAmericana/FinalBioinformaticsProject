#################MYTEST###############################

import os
import subprocess
from scripts.getConcatAlignNamedFasta import  fastas_Aligner, fastas_Concatenator


def test_Aligned_Fastas():
    fastas_Aligner("TestSuit/named_Fastas", "TestSuit/alignedFastas")
    assert os.path.exists("TestSuit/alignedFastas")


def test_Concatenator():
    fastas_Concatenator('TestSuit/alignedFastas')
    assert os.path.isfile("concat.fasta")


