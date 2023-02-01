from scripts.getML_Tree import getBestModel, raxMLExecutor
import os

def test_getBest():
    assert os.path.isfile("TestSuit/model.txt")

def test_raxMLEx():
    lastLine = "23  +I"
    raxMLExecutor(lastLine)
    assert os.path.isfile("RAxML_bootstrap.nwk")