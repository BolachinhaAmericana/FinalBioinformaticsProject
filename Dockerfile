#Made by Silva

FROM snakemake/snakemake

RUN apt-get update
RUN apt install software-properties-common --yes
RUN apt-get install --yes build-essential
RUN pip3 install requests biopython


WORKDIR /lab
COPY Snakefile /lab
COPY scripts /lab
