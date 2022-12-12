#get eSearch result as Argument (dictionary)



eSearchResult= {'Count': '1', 'RetMax': '1', 'RetStart': '0', 'QueryKey': '1', 'WebEnv': 'MCID_63965cfc076e8b4fa25f4ff0', 'IdList': ['9606'], 'TranslationSet': [], 'TranslationStack': [{'Term': 'homo sapiens[All Names]', 'Field': 'All Names', 'Count': '1', 'Explode': 'N'}, 'GROUP'], 'QueryTranslation': 'homo sapiens[All Names]'}

webEnv = eSearchResult["WebEnv"]
queryKey = eSearchResult["QueryKey"]

print(webEnv, queryKey)
 

