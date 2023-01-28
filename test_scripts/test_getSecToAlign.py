import pytest
from scripts.getSecToAlign import createSquences_Fasta_Dir, getFastas


def test_createSquences_Fasta_Dir():
    # testa se a função cria corretamente o diretório
    assert createSquences_Fasta_Dir('Squences_Fasta') == 'Squences_Fasta'
    # testa se a função exclui corretamente o diretório caso já exista
    assert createSquences_Fasta_Dir('Squences_Fasta') == None

def test_getFastas():
    # testa se a função retorna o caminho correto para o diretório
    assert getFastas('Squences_Fasta', db='Nucleotide') == 'Squences_Fasta'
    # testa se a função escreve corretamente os arquivos fasta
    assert getFastas('Squences_Fasta', db='Nucleotide') == None