#!/bin/bash

mkdir ./Downloading...
cd ./Downloading...

#Download and extract taxonkit
wget https://github.com/shenwei356/taxonkit/releases/download/v0.14.0/taxonkit_linux_amd64.tar.gz
tar -zxvf *.tar.gz

#install taxonkit without root 
mkdir -p $HOME/bin/; cp taxonkit $HOME/bin/

#Checking if taxonkit is installed
taxonkit version
taxonkit genautocomplete

#Download and extract database
wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

#setting up env and installing database
mkdir -p $HOME/.taxonkit
cp names.dmp nodes.dmp merged.dmp delnodes.dmp $HOME/.taxonkit

#deleting junk files
cd ..
rm -r ./Downloading...
