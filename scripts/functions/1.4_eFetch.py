# Takes a database and querykey/webenv as args and returns the result of the search as a .xml string
from Bio import Entrez

def eFetch(db, webEnv, queryKey):
    '''
    fetches data using Entrez API

    Arguments:
        takes one NCBI database, the querykey and the web environment from a eSearch
    Vars -
        fetchHandle: Uses Entrez API to fetch information with specified paramethers
        eFetchResult: Reads previous value
    Returns:
        eFetchResult
    '''
    if db== 'Taxonomy':
        fetchHandle = Entrez.efetch(db=db, webenv=webEnv, query_key=queryKey)
        eFetchResult = fetchHandle.read()
    return eFetchResult

'''
Output example:
A unfilthered file with a ton of characters.
'''