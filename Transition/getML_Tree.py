'''get maximul likelyhood tree'''
#!/usr/bin/env python3
import subprocess

def get_best_model(fasta):
    '''finds the best model'''
    subprocess.run(f'modeltest-ng -d nt -i {fasta} -o model -p 4 -f e -h i -s 11 > model.txt', shell = True)

    with open("model.txt", encoding='utf8') as doc:
        content = doc.readlines()
        last_line = content[-5]
    return last_line

def raxml_executor(last_line):
    '''executes raxml'''
    if "+I" in last_line:
        subprocess.run("raxmlHPC -m GTRCATI -p 1234567 -x 1234567 -# autoFC -s concat.fasta -n nwk", shell = True)
    else:
        subprocess.run("raxmlHPC -m GTRCAT -p 1234567 -x 1234567 -# autoFC -s concat.fasta -n nwk", shell = True)
    #subprocess.run("rm -r model*", shell = True)
    #subprocess.run("RAxML_info.nwk", shell = True)

if __name__ == "__main__":
    FINAL_LINE = ""
    get_best_model("concat.fasta")
    raxml_executor(FINAL_LINE)
