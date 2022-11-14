#!/usr/bin/env python3

import subprocess
import sys

ficheiro=sys.argv[1]


subprocess.Popen(f"./fasta_to_Nexus.py {ficheiro}", shell=True)


