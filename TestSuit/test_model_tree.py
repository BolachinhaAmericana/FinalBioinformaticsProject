from scripts.getML_Tree import raxml_executor
import os

def test_getBest():
    assert os.path.isfile("TestSuit/model.txt")

def test_raxMLEx():
    lastLine = "23  +I"
    raxml_executor(lastLine)
    assert os.path.isfile("RAxML_bootstrap.nwk")