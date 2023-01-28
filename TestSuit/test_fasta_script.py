#################MYTEST###############################


from scripts.readconcandalignedfastas import  fastas_Aligner, getNamed_fastas, fastas_Concatenator

def test_Aligned_Fastas():
    fastas_Aligner("named_Fastas", "alignedFastas")

def test_getNamed_Fastas():
    getNamed_fastas("named_Fastas")

def test_Concatenator():
    fastas_Concatenator('alignedFastas')


