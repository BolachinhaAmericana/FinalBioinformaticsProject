

# Accepts output from eSearch function and gets the TaxonId from that query

def getTaxonId(eSearchResult):
    '''
    Will get the Taxon Id of our input search term

    Arguments: 
        Takes eSearch result as argument.
    Vars - 
        taxonId: filthering string
    Returns:
        Filthered Taxon Id
    '''
    taxonId = str(eSearchResult["IdList"])
    taxonId = taxonId.replace("[\'",'')
    taxonId = taxonId.replace("\']",'')
    return taxonId

eSearchResult= {'Count': '1', 'RetMax': '1', 'RetStart': '0', 'QueryKey': '1', 'WebEnv': 'MCID_63965cfc076e8b4fa25f4ff0', 'IdList': ['9606'], 'TranslationSet': [], 'TranslationStack': [{'Term': 'homo sapiens[All Names]', 'Field': 'All Names', 'Count': '1', 'Explode': 'N'}, 'GROUP'], 'QueryTranslation': 'homo sapiens[All Names]'}

taxonId= getTaxonId(eSearchResult)
print(taxonId)

'''
Example Output:
9606
'''