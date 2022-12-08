from Bio import Entrez
import os
os.system(f"rm -r yamete.txt")
os.system(f"touch yamete.txt") 
Entrez.email = "guilherme.vcc.silva@gmail.com"  
search_handle = Entrez.esearch(db="Taxonomy", term="primates[subtree]", usehistory="y", idtype="acc")
search_results = Entrez.read(search_handle)
search_handle.close()

acc_list = search_results["IdList"]
print(len(acc_list))
count = int(search_results["Count"])

webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

batch_size = 9999

for start in range(0, count, batch_size):
    end = count
    print("Going to download record %i to %i" % (start + 1, end))
    fetch_handle = Entrez.efetch(
        db="Taxonomy",
        rettype="gb",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=webenv,
        query_key=query_key,
        idtype="acc",
        )
    data = fetch_handle.read()
    print(data)
    #os.system(f"echo {data} >> yamete.txt")    
    fetch_handle.close()