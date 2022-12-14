FROM snakemake/snakemake

RUN apt-get update
RUN apt-get install --yes build-essential

RUN pip3 install requests biopython
RUN pip3 install taxoniq


WORKDIR /lab
COPY Snakefile /lab
COPY scripts /lab
    

