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

if command -v docker >/dev/null 2>&1; then
    docker pull peixeaquatico/magic_phylogenetic:latest && 
    docker run -it --rm -v $(pwd):/lab -e term="'$term'" -e rank="$rank" -e proximity="$proximity" -e similarity="$similarity" peixeaquatico/magic_phylogenetic:latest bash -c "snakemake --cores all; exec /bin/bash"
else
    echo "Docker is not installed Check this link https://docs.docker.com/desktop/install/linux-install/"
fi
