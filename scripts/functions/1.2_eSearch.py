from Bio import Entrez

# Faz uma query no NCBI

# Aceita uma bd Manual e um term input


def eSearch(db, term):
    Entrez.email = 'noWarningPlease@gmail.com'
    '''
    Performs a search query on NCBI 

    Arguments: 
        Takes the NCBI database type and search term as arguments.
    Vars -
        eSearch: Uses the Entrez API to search our query with history and a retmax of 9999
        eSearchResult: reads eSearch in a variable.
    Returns:
        eSearchResult
    '''
    eSearch = Entrez.esearch(db=db, term=term, usehistory="y", retmax = 9999)
    eSearchResult = Entrez.read(eSearch)
    return eSearchResult

eSearchResult= eSearch('Taxonomy', 'Homo sapiens')
print(eSearchResult)

'''
Example Output.

{'Count': '1', 'RetMax': '1', 'RetStart': '0', 'QueryKey': '1', 'WebEnv': 'MCID_63965cfc076e8b4fa25f4ff0', 'IdList': ['9606'], 'TranslationSet': [], 'TranslationStack': [{'Term': 'homo sapiens[All Names]', 'Field': 'All Names', 'Count': '1', 'Explode': 'N'}, 'GROUP'], 'QueryTranslation': 'homo sapiens[All Names]'}
'''