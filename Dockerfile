FROM snakemake/snakemake

RUN apt-get update
RUN apt-get install --yes build-essential

RUN pip3 install requests biopython
RUN pip3 install taxoniq

RUN wget https://github.com/shenwei356/taxonkit/releases/download/v0.14.0/taxonkit_linux_amd64.tar.gz
RUN tar -zxvf *.tar.gz 
RUN cp taxonkit /usr/local/bin/ 
RUN mkdir -p $HOME/bin/; cp taxonkit $HOME/bin/

#RUN wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
#RUN tar -zxvf taxdump.tar.gz
#RUN mv *.dmp /root/.taxonkit


WORKDIR /lab
COPY Snakefile /lab
COPY scripts /lab
    