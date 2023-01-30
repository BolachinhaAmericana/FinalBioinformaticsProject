#Made by Silva
FROM snakemake/snakemake

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
  
RUN apt-get install software-properties-common --yes
RUN pip3 install requests biopython
RUN apt-get install mafft --yes


WORKDIR /root
RUN git clone https://github.com/ddarriba/modeltest
RUN mkdir build && cd build 
RUN cmake .. && make 
RUN export PATH=$PATH:/modeltest/bin/modeltest-ng /bin

RUN apt install seqmagick --yes
#RUN conda config --add channels bioconda
RUN pip install biopython==1.77
#RUN conda install -c bioconda modeltest-ng
RUN apt-get install mrbayes


WORKDIR /lab
COPY Snakefile /lab
COPY scripts /lab