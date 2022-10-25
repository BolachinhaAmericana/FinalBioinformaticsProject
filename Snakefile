
rule getSecBdTerm:
    input:
        "input.txt"       
    output:
        "file.txt"
    shell:
        'python3 scripts/getSecBdTerm.py {input}'


