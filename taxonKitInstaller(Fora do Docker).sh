#!/bin/bash


# This script must be run with root privileges, if the user is not root. It will simply exit and do nothing.
# This if statement is merely checking whether or not the script is running with root privileges.
if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

mkdir ./download
cd download
mkdir ./1
cd 1
wget https://github.com/shenwei356/taxonkit/releases/download/v0.14.0/taxonkit_linux_amd64.tar.gz
tar -zxvf *.tar.gz 
cp taxonkit /usr/local/bin/
mkdir -p $HOME/bin/; cp taxonkit $HOME/bin/
cd ..
mkdir ./2
cd 2
wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
mv *.dmp /root/.taxonkit
cd ..
cd ..
rm -r download
