FROM python:3.8.2

RUN apt-get update

RUN pip3 install requests biopython
    