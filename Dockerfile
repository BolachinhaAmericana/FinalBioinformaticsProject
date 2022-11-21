FROM snakemake/snakemake

RUN apt-get update

RUN pip3 install requests biopython &&
pip3 install taxoniq

WORKDIR /lab
COPY Snakefile /lab
COPY scripts /lab
    