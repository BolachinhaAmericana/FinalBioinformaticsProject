#!/usr/bin/python
import toytree # Uma biblioteca para o desenho gráfico de árvores           
import toyplot.pdf
import sys

def get_user_args(specieName):
    specieName = sys.argv[1]
    specieName = specieName.replace(" ", "_")
    print(specieName) 
    return specieName


def oppen_tree_file(treeFile):
    with open(treeFile,'r') as file:
        treeFile = file.read()
    tre = toytree.tree(treeFile) 
    return tre

def style_MLtree(specieName,tre):
     # Uma lista com as diversas cores para os nomes das espécies da árvores
    colorlist = ["darkcyan" if specieName in t else "darkorange" for t in tre.get_tip_labels()]
    # Define os estilos a utilizar
    mystyle = {
        "edge_type": 'p',
        "edge_style": {
            "stroke": toytree.colors[2],
            "stroke": "darkcyan",
            "stroke": "darkorange",
            "stroke-width": 2.5,
        },
        "tip_labels_align": True,
        "tip_labels_colors" : colorlist,
        "tip_labels_style": {
            "font-size": "10px",
        },
        "node_labels" : tre.get_node_values("support",1,0),
        "node_sizes" : 7,
        "node_colors": toytree.colors[2],
    }
    # Guarda os valores que retornam do desenho da arvore em canvas, axes, mark
    mlCanvas = tre.draw(height=900,**mystyle)
    return mlCanvas

def style_MBtree(specieName,tre):
    # Uma lista com as diversas cores para os nomes das espécies da árvores
    colorlist = ["darkcyan" if specieName in t else "indigo" for t in tre.get_tip_labels()]


    # Define os estilos a utilizar
    mystyle = {
        "edge_type": 'p',
        "edge_style": {
            "stroke": toytree.colors[2],
            "stroke-width": 2.5,
        },
        "tip_labels_align": True,
        "tip_labels_colors": colorlist,
        "tip_labels_style": {
            "font-size": "10px"
        },
        "node_labels" : tre.get_node_values("support",1,0),
        "node_sizes" : 7,
        "node_colors": toytree.colors[2],
    }


    # Guarda os valores que retornam do desenho da arvore em canvas, axes, mark
    mbCanvas = tre.draw(height=900,**mystyle)
    return mbCanvas

def write_to_pdf(mlCanvas,mbCanvas):
    toyplot.pdf.render(mlCanvas, "tree-plot_bootstrap_ML.pdf")
    toyplot.pdf.render(mbCanvas, "tree-plot_MB.pdf")

if __name__ == '__main__':
    specieName = get_user_args(specieName='')
    tre = oppen_tree_file(treeFile="RAxML_bootstrap.nhk")
    mlCanvas = style_MLtree(specieName,tre)
    tre = oppen_tree_file(treeFile="MyBayes.tre")
    mbCanvas = style_MBtree(specieName,tre)
    write_to_pdf(mlCanvas,mbCanvas)





