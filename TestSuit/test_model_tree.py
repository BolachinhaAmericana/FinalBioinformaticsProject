from scripts.getML_Tree import getBestModel, raxMLExecutor
import os

def test_getBest():
    getBestModel("concat.fasta")
    assert os.path.isfile("model.txt")

def test_raxMLEx():
    lastLine = "23  +I"
    raxMLExecutor("teste_ali.fasta",lastLine)
    assert os.path.isfile("RAxML_bootstrap.nwk")