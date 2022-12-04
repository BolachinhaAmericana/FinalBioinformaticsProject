import taxoniq

file = open('ListTaxonIds.txt')
#Lines = file.readlines

lista = []

for row in file:
    lista.append(row)

goodList = [x[:-1] for x in lista]

for row in goodList:
    t = taxoniq.Taxon(row)
    print(t.scientific_name)
