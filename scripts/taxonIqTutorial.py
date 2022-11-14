import taxoniq

t=taxoniq.Taxon(9606)
nomeEspecie = t.scientific_name
print(f"Scientific Name: {nomeEspecie}")
nomeComum = t.common_name
print(f"Common Name: {nomeComum}")

listaRanks = t.ranked_lineage
#print(f"Ranked Lineage: {listaRanks}")
linhagemCompleta = t.lineage
len(linhagemCompleta) # nº ancestrais conhecidos

relacaoDeRankOrganismo = [(t.rank.name, t.scientific_name) for t in listaRanks]
#print(relacaoDeRankOrganismo)

subespeciesOrganismo = [(c.rank.name, c.common_name) for c in t.child_nodes]
#print(subespeciesOrganismo)

listaAccnRepresentantesGenoma = t.refseq_representative_genome_accessions[:10]
print(listaAccnRepresentantesGenoma)