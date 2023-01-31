FROM debian:buster-slim AS builder
WORKDIR /modeltest

RUN apt update \
    && apt install -y \
    autoconf \
    automake \
    libtool \
    bison \
    flex \
    cmake \
    build-essential \
    git
RUN git clone https://github.com/ddarriba/modeltest 
RUN cd modeltest && mkdir build && cd build \
    && cmake .. && make


FROM snakemake/snakemake

RUN apt update
RUN apt-get install software-properties-common --yes
RUN python3 -m pip install requests biopython
RUN apt-get install mafft --yes

COPY --from=builder /modeltest/modeltest/bin/modeltest-ng /bin

RUN python3 -m pip install seqmagick
RUN apt-get install mrbayes --yes
RUN python3 -m pip install toytree
RUN apt-get install raxml


WORKDIR /lab
COPY Snakefile /lab
COPY scripts /lab
