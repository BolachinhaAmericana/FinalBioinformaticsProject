#!/bin/bash

while getopts "s:r:p:i:" opt; do
    case ${opt} in
        s ) term="$OPTARG" ;;
        r ) rank="$OPTARG" ;;
        p ) proximity="$OPTARG" ;;
        i ) similarity="$OPTARG" ;;
        \? ) echo "Error: Invalid option -$OPTARG"
             echo "Usage: executer.sh -s <term> -r <rank> -p <proximity> -i <similarity>"
             exit 1
    esac
done

if [ -z "$term" ] || [ -z "$rank" ] || [ -z "$proximity" ] || [ -z "$similarity" ]; then
    echo "Error: Please provide all required arguments" 
    echo "Usage: executer.sh -s <term> -r <rank> -p <proximity> -i <similarity>"
    exit 1
fi

#docker build -t myimage . && 
#docker run --rm -v $(pwd):/data myimage 
term="'$term'" rank="$rank" proximity="$proximity" similarity="$similarity"  snakemake --cores all
