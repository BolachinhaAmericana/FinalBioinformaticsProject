import sys
import taxoniq


#def obterArgumentosDoUtilizador():
taxonomia = ""
rankTaxonomia = sys.argv[1]
term = sys.argv[2]
nome = taxoniq.Taxon(scientific_name=term)
#return rankTaxonomia, term



#nomeEspecie = t.scientific_name
#print(f"Scientific Name: {nomeEspecie}")
#nomeComum = t.common_name
#print(f"Common Name: {nomeComum}")

#linhagemCompleta = t.lineage
#print(f"nº ancestrais conhecidos: {len(linhagemCompleta)}") # nº ancestrais conhecidos

listaRanks = nome.ranked_lineage
#print(f"Ranked Lineage: {listaRanks}")


relacaoDeRankOrganismo = [(nome.rank.name, nome) for nome in listaRanks]
print(f"Relacao De Rank Organismo: {relacaoDeRankOrganismo}")


for i in relacaoDeRankOrganismo:
    if rankTaxonomia in str(i):
        taxonomia = i[1]
        print(type(taxonomia))
        print(taxonomia.scientific_name)
        print(taxonomia.common_name)
        print([(taxonomia.rank.name, taxonomia.scientific_name) for taxonomia in taxonomia.ranked_lineage])
        print(taxonomia.url)
        break

#metes o taxonomia.scientific_name do for como argumento term no getSecBdTerm
#da me ai o check valente para depois meter no snakemake + docker

#Com lepra a baixo
#listaAccnRepresentantesGenoma = taxonomia.refseq_representative_genome_accessions[:10]
#print(f"ACCNs: {listaAccnRepresentantesGenoma}")
