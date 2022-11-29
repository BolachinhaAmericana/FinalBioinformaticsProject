FROM snakemake/snakemake

RUN apt-get update
RUN apt-get install --yes build-essential
<<<<<<< HEAD

RUN pip3 install requests biopython
RUN pip3 install taxoniq

RUN wget https://github.com/shenwei356/taxonkit/releases/download/v0.14.0/taxonkit_linux_amd64.tar.gz
RUN tar -zxvf *.tar.gz 
RUN cp taxonkit /usr/local/bin/ 
RUN mkdir -p $HOME/bin/; cp taxonkit $HOME/bin/

#RUN wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
#RUN tar -zxvf taxdump.tar.gz
#RUN mv *.dmp /root/.taxonkit

=======
RUN pip3 install requests biopython
RUN pip3 install taxoniq
>>>>>>> 93aab84c1968bfd0c3d0a7d4725907c00a90c4e1

WORKDIR /lab
COPY Snakefile /lab
COPY scripts /lab
    

