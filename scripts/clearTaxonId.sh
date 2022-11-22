#!/bin/bash

sed 's/\[//' ./1.txt > ./2.txt
sed 's/\]//' ./2.txt > ./3.txt

cat ./3.txt

rm ./1.txt
rm ./2.txt
rm ./3.txt
