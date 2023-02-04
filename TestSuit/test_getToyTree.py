from scripts.get_toy_tree import style_ml_tree, mb_tree_styler, write_to_pdf
import toytree
import toyplot.pdf
import os

def test_style_ml_tree():
    specie_name = "A"
    tree = toytree.tree("(A:1,B:2);")
    result = style_ml_tree(specie_name, tree)
    assert isinstance(result, toyplot.canvas.Canvas)

def test_mb_tree_styler():
    specie_name = "A"
    tree = toytree.tree("(A:1,B:2);")
    result = mb_tree_styler(specie_name, tree)
    assert isinstance(result, toyplot.canvas.Canvas)

def test_write_to_pdf():
    ml_canvas = toyplot.canvas.Canvas()
    mb_canvas = toyplot.canvas.Canvas()
    write_to_pdf(ml_canvas, mb_canvas)
    # Verifica se os arquivos PDF foram criados corretamente
    assert os.path.isfile("tree-plot_bootstrap_ML.pdf")
    assert os.path.isfile("tree-plot_MB.pdf")