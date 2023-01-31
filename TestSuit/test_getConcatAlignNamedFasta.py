#################MYTEST###############################

import os
import shutil
from scripts.getConcatAlignNamedFasta import  fastas_Aligner, getNamed_fastas, fastas_Concatenator

def test_Aligned_Fastas():
    fastas_Aligner("named_Fastas", "alignedFastas")
    assert os.path.exists("alignedFastas")

def test_getNamed_Fastas():
    dir = "named_Fastas"
    getNamed_fastas(dir)
    assert os.path.exists(dir)
    assert len(os.listdir(dir)) > 0

def test_Concatenator():
    fastas_Concatenator('alignedFastas')
    assert os.path.isfile("concat.fasta")
    


