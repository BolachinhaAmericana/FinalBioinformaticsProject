#!/usr/bin/env python

import subprocess
import sys

ficheiro=sys.argv[1]


subprocess.Popen(f"./fasta_to_Nexus.py {ficheiro}", shell=True)


