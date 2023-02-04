#!/usr/bin/python
'''Tree design'''

import sys
import toytree # Uma biblioteca para o desenho gráfico de árvores
import toyplot.pdf
from essentials import get_user_arguments

def open_tree_file(tree_file):
    '''Opens tree file'''
    with open(tree_file,'r', encoding='utf8') as text:
        tree_file = text.read()
    tre = toytree.tree(tree_file)
    print(tre)
    return tre

def style_ml_tree(specie_name,tre):
    '''sets tree style'''
     # Uma lista com as diversas cores para os nomes das espécies da árvores
    color_list = ["darkcyan" if specie_name in t else "darkorange" for t in tre.get_tip_labels()]
    # Define os estilos a utilizar
    my_style = {
        "edge_type": 'p',
        "edge_style": {
            "stroke": toytree.colors[2],
            "stroke": "darkcyan",
            "stroke": "darkorange",
            "stroke-width": 2.5,
        },
        "tip_labels_align": True,
        "tip_labels_colors" : color_list,
        "tip_labels_style": {
            "font-size": "10px",
        },
        "node_labels" : tre.get_node_values("support",1,0),
        "node_sizes" : 7,
        "node_colors": toytree.colors[2],
    }
    # Guarda os valores que retornam do desenho da arvore em canvas, axes, mark
    ml_canvas ,axes, mark = tre.draw(height=900,**my_style)
    return ml_canvas
def mb_tree_styler(specie_name, tre):
    ''' sets the style for the mb tree'''
    # Uma lista com as diversas cores para os nomes das espécies da árvores
    color_list = ["darkcyan" if specie_name in t else "indigo" for t in tre.get_tip_labels()]
    # Define os estilos a utilizar
    my_style = {
        "edge_type": 'p',
        "edge_style": {
            "stroke": toytree.colors[2],
            "stroke-width": 2.5,
        },
        "tip_labels_align": True,
        "tip_labels_colors": color_list,
        "tip_labels_style": {
            "font-size": "10px"
        },
        "node_labels" : tre.get_node_values("support",1,0),
        "node_sizes" : 7,
        "node_colors": toytree.colors[2],
    }

    # Guarda os valores que retornam do desenho da arvore em canvas, axes, mark
    mb_canvas, axes, mark = tre.draw(height=900,**my_style)
    return mb_canvas

def write_to_pdf(ml_canvas,mb_canvas):
    '''writes input to pdf file'''
    toyplot.pdf.render(ml_canvas, "tree-plot_bootstrap_ML.pdf")
    toyplot.pdf.render(mb_canvas , "tree-plot_MB.pdf")

if __name__ == '__main__':
    ml_tree, mr_bayes, specie = get_user_arguments(3)
    specie = specie.replace(" ", "_")

    tree = open_tree_file(tree_file=ml_tree)
    max_likelyhood_canvas = style_ml_tree(specie,tree)
    tree = open_tree_file(tree_file=mr_bayes)
    mr_bayes_canvas = mb_tree_styler(specie,tree)
    write_to_pdf(max_likelyhood_canvas,mr_bayes_canvas)
