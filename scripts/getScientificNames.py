import taxoniq


#Something is wrong here!!!


file = open('ListTaxonIds.txt')
#Lines = file.readlines

for row in file:
    t = taxoniq.Taxon(row)
    print(t.scientific_name)