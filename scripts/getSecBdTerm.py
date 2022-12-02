#! /usr/bin/python
import sys
from Bio import Entrez


def getArgumentFromUser():
    ''' 

    Obtem os argumentos "Data Base" e o "term" escritos pelo utilizador
    return: 
    Os argumentos "Data Base" db (String) e o "term" term (String) escritos pelo utilizador

    '''

    db = sys.argv[1]
    term = sys.argv[2]
    return db, term


def getResultsFromESearch(db, term):
    ''' 

    Obtem o resultado de um ESearch feito à "Entrez API" com a "Data Base" e o "term" escritos pelo utilizador
    param: 
    db ,term: "Data Base" e o "term" escritos pelo utilizador
    return: 
    o resultado da pesquisa feito ESearch feito à "Entrez API"

    '''

    eSearch = Entrez.esearch(db=db, term=term, usehistory="y")
    result = Entrez.read(eSearch)
    return result


def getWebEnv(result):
    ''' 

    Obtem o WebEnv a partir do resultado da função anterior
    param: 
    resultado: o resultado da pesquisa feito ESearch feito à "Entrez API"
    return: 
    O WebEnv do resultado

    '''

    webEnv = result["WebEnv"]
    return webEnv


def getQueryKey(result):
    ''' 

    Obtem a QueryKey a partir do resultado da primeira função
    param: 
    resultado: o resultado da pesquisa feito ESearch feito à "Entrez API"
    return: 
    A QueryKey do resultado

    '''

    queryKey = result["QueryKey"]
    return queryKey


def getInformactionEFetch(queryKey, webEnv):
    ''' 

    Obtem a informação do EFetch a partir da "Data Base", da QueryKey e do WebEnv que veem das funções anteriores
    param: 
    queryKey, webEnv: a QueryKey e o WebEnv do resultado feito com o ESearch feito à "Entrez API"
    return: 
    A informacao pelo feito EFetch feiro com a "Data Base", QueryKey e WebEnv no formato fasta

    '''

    fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey, rettype='xml')
    informaction = fetchHandle.read()
    return informaction

if __name__ == '__main__':
    db, term = getArgumentFromUser()
    result = getResultsFromESearch(db, term)
    queryKey, webEnv = getQueryKey(result), getQueryKey(result)
    informaction = getInformactionEFetch(queryKey, webEnv)
    print(informaction)
    sys.stderr.write("Foi adicionado com sucesso um ficheiro sequenciacao.Fasta com o output >> ")
    
